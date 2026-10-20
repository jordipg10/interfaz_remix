!> \file aqueous_chemistry_m.f90
!> \brief Defines the aqueous_chemistry_c type and its associated procedures for aqueous phase chemical modeling.
!>
!> This module provides the core data structures and algorithms for simulating aqueous chemistry, including species concentrations, activities, activity coefficients, and coupling with solid and gas phases. It supports equilibrium and kinetic reactions, mixing iterations, and advanced speciation algorithms for reactive transport.
!>
!> The reactive-mixing iteration routines come in explicit (Euler) and implicit variants and in
!> lumped and consistent forms. The implicit consistent routines (compute_r_tilde_impl_opt1..4)
!> accept an optional upstream `curr_water`, so a cell can be coupled to its already-updated
!> upstream neighbours; this is what the 2D consistent WMA solver relies on.
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2024
module aqueous_chemistry_m
    use mineral_zone_m, only: mineral_zone_c                              !< Import mineral zone class for mineral interactions
    use aq_phase_m, only: aq_phase_c                             !< Import aqueous phase class for species and reactions
    use local_chemistry_m, only: local_chemistry_c                         !< Import base local chemistry class for inheritance
    use solid_chemistry_m, only: solid_chemistry_c !< Import solid phase chemistry classes for mineral interactions
    use reactive_zone_m, only: reactive_zone_c                              !< Import reactive zone class (zone argument of the within-zone mixing routines)
    use gas_chemistry_m, only: gas_chemistry_c                          !< Import gas chemistry for gas-water equilibrium
    use arrays_m, only: int_array_c, diag_matrix_c, id_matrix, norm_mat_inf
    use metodos_sist_lin_m, only: LU_lin_syst
    use vectors_m, only: inf_norm_vec_real, outer_prod_vec
    use params_aq_sol_m, only: params_aq_sol_s                             !< Import aqueous solution parameter structure
    use array_ops_m, only: get_ind_int, is_int_in_1D_array                 !< Import array utility functions for index operations
    use reaction_m, only: reaction_c                                        !< Import base reaction class (per-reaction T-dependent logK helper)
    implicit none                                                           !< Require explicit variable declaration (Fortran safety)
    save   !< Preserve module variables between procedure calls      
    private !< Make all entities private by default
    public ::  initialise_iterative_method                     !< Expose iterative method initialisation subroutine
    public :: mixing_iter_comp, reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt4, mixing_iter_comp_ideal, compute_r_tilde_impl_opt1, reactive_mixing_iter_EI_eq_anal_ideal, &
        mixing_iter_species, reactive_mixing_iter_EE_kin, reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_c_mix, compute_c_mix_global_cons, compute_c_mix_global_lump, compute_r_tilde_impl_opt2, compute_r_tilde_impl_opt3, &
        compute_r_tilde_impl_opt4_idx, compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, &
        compute_r_tilde_impl_opt3_idx, &
        compute_c_hat_patankar, add_limited_increment, compare_old_reaction_zones, mix_within_zone_new_rz, &
        compute_c_mix_by_zones_cons, compute_c_mix_by_zones_lump, &
        add_curr_step_react_mix_contrib_opt1, add_curr_step_react_mix_contrib_opt2, &
        add_curr_step_react_mix_contrib_opt3, &
        write_celldbg_legend, write_celldbg_cell, reactive_mixing_iter_EI_eq_kin_anal_opt2, &
        reactive_mixing_iter_RKF45_eq_kin_anal_opt2 !< Expose mixing subroutines for external use
    !> \brief Aqueous chemistry class for modeling chemical reactions in the aqueous phase.
    !>
    !> This type extends local_chemistry_c and encapsulates all properties and methods required for aqueous chemical speciation, activity calculations, and coupling with solid and gas chemistry. It manages concentrations, activities, activity coefficients, and provides interfaces for mixing, reaction iteration, and property calculation.
    type, public, extends(local_chemistry_c) :: aqueous_chemistry_c
        type(params_aq_sol_s) :: params_aq_sol                             !< Aqueous solution parameters
        real(kind=8) :: ionic_strength                                     !< Ionic strength (I) = 0.5 * sum(m_i * z_i^2) for activity coefficient calculations
        real(kind=8) :: pH                                             !< pH = -log10[H+] - measure of hydrogen ion activity (acidity); stays 0 when h+ is not in the aqueous phase
        real(kind=8) :: pe                                                 !< pe = -log10[e-] - measure of electron activity (redox potential)
        real(kind=8) :: salinity                                           !< Total dissolved solids concentration (salinity measure)
        real(kind=8) :: alkalinity                                         !< Alkalinity = m_OH- + m_HCO3- + 2*m_CO3-2 (acid neutralizing capacity)
        real(kind=8), allocatable :: conc_mix(:)                          !< Aqueous variable-activity species concentrations of the water's CONSERVATIVE mix of the current step (previous-step concentrations mixed; for the by-zones mixers including the previous-step kinetic term of the within-zone sub-mixtures), ordered [aqueous primary (1:num_aq_prim_species) | aqueous secondary var act (num_aq_prim_species+1:num_aq_var_act_species)]. Stored by the mixers (compute_c_mix_cons_global, cross_zone_mix_core). Written to the 'conc_mix' output blocks.
        real(kind=8), allocatable :: conc_hat(:)                          !< Aqueous variable-activity species concentrations of the mixed waters that ENTER THE NEWTON SPECIATION of the current step, same layout as conc_mix [aq prim | aq sec var act]: c_hat = c_hat_old + Delta_t*r_tilde, i.e. the mix including the previous-step (1-theta_r) reaction terms of all mixing waters and the current-step (theta_r) reaction terms of the mixing waters EXCEPT this water's own, which the speciation solves for (mixing_iter_comp(_ideal): c_hat is its input). Stored by the reactive mixing solvers once c_hat is formed (set_conc_hat); for a within-zone sub-mixture, its c_hat_z (mix_within_zone_core). Written to the 'conc_hat' output blocks.
        real(kind=8), allocatable :: logK(:)                              !< Per-water equilibrium constants K_eq(T) of the reactive zone's equilibrium reactions, evaluated at THIS water's temperature (state variable). The shared reaction objects hold only the reference-T eq_cst; this T-dependent value is built per-water by compute_log_K_aq_chem (via compute_logK_dep_T).
        real(kind=8), allocatable :: logK_tilde(:)                        !< Per-water transformed eq. constants (flag_comp=.false., all secondary species): inv_Se_2Ã‚Â·log10(K), with constant-activity gas partial pressures folded into K (K/P_gas). Kept per-water (NOT on the shared reactive zone) because waters with different partial pressures may share a reactive zone.
        real(kind=8), allocatable :: logK_star(:)                         !< Per-water transformed eq. constants (flag_comp=.true., variable activity secondary species): inv_Se_nc_2Ã‚Â·log10(K)
        class(solid_chemistry_c), pointer :: solid_chemistry=>null()       !< Pointer to associated solid chemistry object (minerals, sorption)
        class(solid_chemistry_c), pointer :: solid_chemistry_old=>null()   !< Pointer to associated solid chemistry object in previous time step (Lagrangian particle tracking)
        class(gas_chemistry_c), pointer :: gas_chemistry=>null()           !< Pointer to associated gas chemistry object (gas-water equilibrium)
        class(aq_phase_c), pointer :: aq_phase=>null()                     !< Pointer to aqueous phase definition (species list, reactions)
        integer(kind=4), allocatable :: indices_aq_species(:)              !< Indices of aqueous species in state variables
        integer(kind=4), allocatable :: indices_aq_phase(:)                !< Indices of aqueous species in the aqueous phase object
        integer(kind=4), allocatable :: ind_prim_species(:)                !< Indices of primary species (basis species for speciation calculations)
        integer(kind=4), allocatable :: ind_prim_species_old(:)                !< Indices of primary species (basis species for speciation calculations)
        integer(kind=4), allocatable :: ind_var_act_species(:)             !< Indices of variable activity species in global state variables (concentrations computed from equilibrium)
        integer(kind=4), allocatable :: ind_sec_species(:)                 !< Indices of secondary species (computed from primary via mass action)
        integer(kind=4), allocatable :: ind_sec_species_old(:)                 !< Indices of secondary species (computed from primary via mass action)
        real(kind=8), allocatable :: saturation_mins(:)                        !< Saturation index Omega = IAP/K(T) [-] of EVERY mineral of the chemical system (kinetic AND equilibrium), indexed by the mineral's position in chem_syst%minerals (kinetic minerals first, so the kinetic index ind_min of compute_saturation_kin_min is also its position here; then the equilibrium minerals). Allocated by allocate_saturation_mins (called on first use) and initialised to 1 (equilibrium); kinetic entries are refreshed by compute_saturation_kin_min at each rate/Jacobian evaluation of the water's current mineral zone, equilibrium entries by compute_saturation_eq_min, and compute_saturation_mins refreshes all of them at the water's current state (end of each time step). Entries never computed keep 1. Per-water state: travels with the parcel in copy_aq_chem_Lagr.
        type(int_array_c) :: indices_rk                                    !< Indices of species relevant for kinetic reaction rate calculations
        type(int_array_c) :: indices_rk_old                                !< Indices of species relevant for kinetic reaction rate calculations in previous time step
    contains                                                                !< Type-bound procedures section begins
        procedure :: read_wat_type_CHEPROO                         !< Read water type from CHEPROO format input files
    !> Set procedures - Methods for setting/updating chemical state
        procedure :: set_concentrations=>set_conc_aq_species       !< Set aqueous species concentrations from external input
        procedure :: set_act_aq_species                            !< Set activities of aqueous species (for equilibrium calculations)
        procedure :: set_act_diss_solids                           !< Set activities of dissolved solid species
        procedure :: set_conc_aq_prim_species                      !< Set concentrations of aqueous primary species only
        procedure :: set_conc_prim_species                         !< Set concentrations of all primary species (aqueous + solid)
        procedure :: set_conc_var_act_species                       !< Set concentrations of variable activity species (equilibrium-dependent)
        procedure :: set_conc_sec_var_act_species                   !< Set concentrations of secondary variable activity species
        procedure :: set_conc_sec_aq_species                        !< Set concentrations of all secondary aqueous species
        procedure :: set_conc_sec_aq_var_act_species                !< Set concentrations of secondary aqueous variable activity species
        procedure :: set_conc_sec_species                           !< Set concentrations of all secondary species (aqueous + solid)
        procedure :: set_indices_aq_species_aq_chem                 !< Set indices mapping aqueous species to global state variables
        procedure :: set_eq_rxn_rates                     !< Set equilibrium reaction rates during a time step
        procedure :: set_spec_alg_aq_chem                           !< Set speciation algebra object for equilibrium calculations
        procedure :: set_ind_species                                !< Set species indices for chemical system navigation
        procedure :: set_ind_species_react_zone                     !< set_ind_species + the reactive zone's primary/secondary swap (zone column order)
        procedure :: set_pH                                         !< Set pH value (negative log of hydrogen ion activity)
        procedure :: set_pe                                         !< Set pe value (negative log of electron activity - redox state)
        procedure :: set_ionic_strength                             !< Set ionic strength for activity coefficient calculations
        procedure :: set_conc_single_species                        !< Set concentration of a single chemical species
        procedure :: set_log_act_coeffs                             !< Set logarithmic activity coefficients for non-ideal solutions
        procedure :: set_solid_chemistry                            !< Associate solid chemistry object (minerals, sorption sites)
        procedure :: set_solid_chemistry_old                        !< Associate old solid chemistry object (minerals, sorption sites)
        procedure :: set_gas_chemistry                              !< Associate gas chemistry object (gas-water equilibrium)
        procedure :: set_aq_phase                                   !< Associate aqueous phase object (species definitions, reactions)
        procedure :: set_indices_rk                                 !< Set indices of species relevant for kinetic reaction calculations
        procedure :: set_Rk                                         !< Set kinetic reaction rates for time integration
        procedure :: set_rxn_rate_kin_avg                                    !< Set mean kinetic reaction rates for time stepping
        procedure :: set_rk_old                                  !< Set kinetic reaction rates from previous time step (for time averaging and predictor-corrector)
        procedure :: set_attributes_dep_rz                     !< Set attributes dependent on reactive zone
        procedure :: set_attributes_dep_mz                     !< Set attributes dependent on mineral zone
    !> Allocate procedures - Memory allocation methods for dynamic arrays
        procedure :: allocate_rxn_rates_and_amounts                        !< Allocate memory for kinetic reaction rate arrays
        procedure :: allocate_conc_aq_species                       !< Allocate memory for aqueous species concentration arrays
        procedure :: allocate_conc_mix                          !< Allocate memory for mixed aqueous species concentration arrays
        procedure :: allocate_conc_hat                          !< Allocate memory for the aqueous species concentrations entering the Newton speciation
        procedure :: set_conc_hat                               !< Store the aqueous blocks of the Newton speciation input c_hat in conc_hat
        procedure :: allocate_contrib_rxn_amnt                      !< Allocate memory for primary species concentration arrays
        procedure :: allocate_activities_aq_species                 !< Allocate memory for aqueous species activity arrays
        procedure :: allocate_log_act_coeffs_aq_chem                !< Allocate memory for logarithmic activity coefficient arrays
        procedure :: allocate_indices_aq_species_aq_chem            !< Allocate memory for aqueous species index arrays
        procedure :: allocate_indices_aq_phase_aq_chem              !< Allocate memory for aqueous phase index arrays
        procedure :: allocate_saturation_mins                       !< Allocate memory for the saturation indices of all minerals of the chemical system (initialised to 1)
    !> Compute procedures - Calculation methods for chemical properties
        procedure :: compute_activities_diss_solids                 !< Compute activities of dissolved solid species
        procedure :: compute_a1                                     !< Compute primary species activities (a1 vector)
        procedure :: compute_a2                                     !< Compute secondary species activities (a2 vector)
        procedure :: compute_activities_aq                          !< Compute all aqueous species activities
        procedure :: compute_activities                             !< Compute activities for all species types
        procedure :: compute_activity                               !< Compute activity of a single species
        procedure :: compute_act_from_MAL                           !< Compute activities from mass action law equilibrium
        procedure :: compute_act_water                              !< Compute water activity (solvent activity)
        procedure :: compute_conc_ideal_water                       !< Compute water concentration assuming ideal solution
        procedure :: compute_molarities                             !< Compute molarities (mol/L) from molalities
        procedure :: compute_molalities                             !< Compute molalities (mol/kg H2O) from concentrations
        procedure :: compute_react_surfaces_per_kg_water            !< Reactive surfaces of the solid chemistry converted from m2 mineral/L water to m2 mineral/kg H2O (molality basis), returned as an argument
        procedure :: compute_react_surfaces_per_L_water             !< Reactive surfaces given in m2 mineral/kg H2O converted to m2 mineral/L water and stored in the solid chemistry
        procedure :: compute_density_water                          !< Compute the pure-water density [kg/L] from the temperature [K] (Kell 1975 correlation) and store it in density
        procedure :: compute_params_aq_sol => compute_params_aq_sol_aq_chem !< Evaluate this water's Debye-Hückel A and B at its temperature and water density
        procedure :: compute_salinity                               !< Compute total dissolved solids (salinity)
        procedure :: compute_alkalinity                             !< Compute alkalinity (acid neutralizing capacity)
        procedure :: compute_pH                                     !< Compute pH from hydrogen ion activity
        procedure :: compute_conc_comp                              !< Compute component concentrations for mass balance
        procedure :: compute_conc_comp_cst_act                      !< Compute component concentrations with constant activity species
        procedure :: compute_re_mean                                    !< Compute mean reaction extent for time averaging
        procedure :: compute_r_tilde_lump                                !< Compute reaction extent with spatial lumping
        procedure :: compute_Re                                     !< Compute reaction extent (equilibrium + theta_r-weighted kinetic term, PDF steps 2d.v-vii)
        procedure :: compute_rxn_rate_kin                                 !< Compute new kinetic reaction rates for current time step
        procedure :: compute_Jacobian_rk_anal                       !< Compute analytical Jacobian of kinetic reaction rates
        procedure :: compute_rk_Jac_rk_anal                         !< Compute kinetic rates and their analytical Jacobian
        procedure :: compute_rk_Jac_rk_incr_coeff                   !< Compute kinetic rates and Jacobian using incremental coefficients
        procedure :: compute_logK_PHREEQC_aq_chem                  !< Compute logarithmic equilibrium constants for aqueous reactions
        procedure :: compute_logK_star      !< Build per-water transformed constant logK_star from the unmodified K
        procedure :: compute_logK_tilde_aq_chem                     !< Build per-water transformed constant logK_tilde from the unmodified K
        procedure :: compute_ionic_strength                         !< Compute ionic strength from species concentrations and charges
        procedure :: compute_d_log_gamma_d_I_aq_chem                !< Compute derivative of log activity coefficients with respect to ionic strength
        procedure :: compute_res_Jac_res_anal                       !< Compute residual and analytical Jacobian for equilibrium solving
        procedure :: compute_res_Jac_res_anal_ideal                 !< Compute residual and Jacobian assuming ideal solution behavior
        procedure :: compute_res_Jac_res_anal_exch                  !< Compute residual and Jacobian for ion exchange reactions
        procedure :: compute_res_Jac_res_incr_coeff                 !< Compute residual and Jacobian using incremental coefficient method
        procedure :: compute_res_init                               !< Compute initial residual for equilibrium iteration
        procedure :: compute_log_act_coeff_wat                     !< Compute logarithmic activity coefficient of water
        procedure :: compute_log_act_coeffs_aq_chem                 !< Compute logarithmic activity coefficient of water
        procedure :: compute_log_act_coeffs_var_act_ideal           !< Compute logarithmic activity coefficient of water
        procedure :: compute_log_act_coeffs_ideal                    !< Compute logarithmic activity coefficients assuming ideal solution behavior
        procedure :: compute_log_Jacobian_act_coeffs_aq_chem        !< Compute logarithmic activity coefficient of water
        procedure :: compute_saturation_kin_min                     !< Compute (and store) the saturation index of a kinetic mineral
        procedure :: compute_saturation_eq_min                      !< Compute (and store) the saturation index of an equilibrium mineral
        procedure :: compute_saturation_mins                        !< Compute (and store) the saturation indices of all minerals of the chemical system
        procedure :: compute_partial_pressures_eq                   !< Compute equilibrium gas partial pressures from the mass-action law (aqueous activities)
        procedure :: compute_rxn_amnt_kin                            !< Compute kinetic reaction rates matrix
        procedure :: compute_contrib_rxn_amnt_kin                        !< Compute contributions to kinetic reaction amounts
        procedure :: compute_contrib_rxn_amnt_eq                         !< Compute contributions to equilibrium reaction amounts (rho_eq = lambda_r*Se_nc^T*Re)
        procedure :: compute_rxn_rate_kin_avg                        !< Compute average kinetic reaction rates during a time step
        procedure :: sum_eq_rxn_rates
    !> Get procedures - Accessor methods for retrieving chemical data
        procedure :: get_rk_new                                     !< Get new kinetic reaction rates (current time step)
        procedure :: get_rk_old                                     !< Get previous-step kinetic reaction rates in the OLD kinetic set (indices_rk_old / solid_chemistry_old; pairs with get_Sk_nc_old)
        procedure :: get_rk_old_new_set                             !< Get previous-step kinetic reaction rates in the NEW kinetic set (indices_rk / current solid_chemistry; pairs with get_Sk_nc)
        procedure :: get_rxn_amnt_kin                                !< Get kinetic reaction amounts
        procedure :: get_rxn_amnt_eq                                 !< Get equilibrium reaction amounts in the reactive zone's stoichiometric-matrix row order
        procedure :: get_Sk_nc                                      !< Get kinetic stoichiometric matrix for variable activity species
        procedure :: get_S_nc                                      !< Get stoichiometric matrix for variable activity species
        procedure :: get_S_nc_old
        procedure :: get_S_nc_non_aq                                !< Get NON-aqueous stoichiometric matrix (current reactive zone) for the Lagrangian previous-step split
        procedure :: get_S_nc_old_aq                                !< Get old aqueous stoichiometric matrix
        procedure :: get_r_old
        procedure :: get_r_old_aq                                   !< Get AQUEOUS-only old reaction rates (matches get_S_nc_old_aq rows)
        procedure :: get_r_old_non_aq                               !< Get NON-aqueous old reaction rates (matches get_S_nc_non_aq rows)
        procedure :: get_r_old_old                                  !< Get full reaction rates from two steps ago (equilibrium rxn_rate_eq_old + kinetic rk_old_old)
        procedure :: get_r_old_old_old                              !< Get full reaction rates from three steps ago (equilibrium rxn_rate_eq_old_old_old + kinetic rk_old_old_old)
        procedure :: get_r_new                                      !< Get kinetic stoichiometric matrix for variable activity species
        procedure :: get_Sk_nc_old                                  !< Get old kinetic stoichiometric matrix for variable activity species
        procedure :: get_c1                                         !< Get primary species concentrations vector
        procedure :: get_conc_old_old_prim_new_rz                                 !< Get primary aqueous species concentrations
        procedure :: get_a1                                         !< Get primary species activities vector
        procedure :: get_conc_var_act_aq_perm                            !< Get non-constant activity aqueous species concentrations (permeable)
        procedure :: get_log_gamma_var_act                          !< Get logarithmic activity coefficients for variable activity species
        procedure :: get_c2                                         !< Get secondary species concentrations vector
        procedure :: get_c2v                                        !< Get secondary non-constant activity species concentrations
        procedure :: get_conc_var_act                                    !< Get all non-constant activity species concentrations
        procedure :: get_conc_old_var_act_new_rz                                !< Get previous-step variable activity concentrations, laid out in the CURRENT reactive zone
        procedure :: get_react_mix_mass_bal_terms                   !< Per-species volume-weighted reactive-mixing mass-balance terms (LHS, transport, reaction)
        procedure :: get_conc                                       !< Get all species concentrations (primary + secondary)
        procedure :: get_eq_csts_aq_chem                            !< Get modified equilibrium constants (constant-activity minerals & gases folded in)
        procedure :: get_log_gamma2v                               !< Get logarithmic activity coefficients for secondary variable activity species
        procedure :: get_log_gamma                                 !< Get logarithmic activity coefficients for all species
        procedure :: get_log_gamma2                                 !< Get logarithmic activity coefficients for all secondary species
        procedure :: get_log_gamma1                                !< Get logarithmic activity coefficients for all primary species
        procedure :: get_u_aq                                       !< Get aqueous component concentrations vector
        procedure :: get_re_mean
    !> Speciation procedures - Chemical equilibrium calculation methods
        procedure :: compute_c2_from_c1_ideal                       !< Compute secondary concentrations from primary activities (ideal)
        procedure :: compute_c2v_from_c1_ideal                     !< Compute secondary non-constant from primary (ideal)
        procedure :: compute_c2v_from_c1_Picard                    !< Compute secondary non-constant from primary (Picard iteration)
        procedure :: compute_c2_from_c1_Picard                      !< Compute secondary from primary (Picard iteration)
        procedure :: compute_c_nc_from_u_Newton                     !< Compute non-constant concentrations from components (Newton method)
        procedure :: compute_c_nc_from_u_Newton_ideal               !< Compute non-constant concentrations (Newton, ideal solution)
        procedure :: compute_res_spec                               !< Compute residual function for equilibrium chemistry
        procedure :: compute_res_spec_cst_act                       !< Compute residual with constant activity coefficients
        procedure :: compute_dc2_dc1_ideal                          !< Jacobian computation for ideal solution assumption
        procedure :: compute_dc2v_dc1_ideal                        !< Non-constant Jacobian for ideal solutions
        procedure :: compute_dc2v_dc1                              !< General non-constant Jacobian computation
        procedure :: compute_dlogc2v_dlogc1                        !< Log-log Jacobian d(log10 c2v)/d(log10 c1), non-ideal (log-space Newton assembly)
        procedure :: compute_dc2_dc1                                !< General Jacobian matrix computation
        procedure :: compute_c_from_u_bin_syst                      !< Concentration from binary system molality
    !> Mixing procedures for water-rock interaction and transport coupling
        procedure :: mixing_iter_comp                               !< Iterative component mixing with activity corrections
        procedure :: mixing_iter_comp_ideal                         !< Component mixing assuming ideal solutions
        procedure :: mixing_iter_species                            !< Species-level iterative mixing procedure
        procedure :: reactive_mixing_iter_EE_eq_kin                    !< Water mixing with explicit equilibrium & kinetics
        procedure :: reactive_mixing_iter_EE_eq_kin_ideal              !< Water mixing EE method (ideal solution)
        procedure :: reactive_mixing_iter_EE_kin               !< EE kinetics-only water mixing (ideal)
        procedure :: reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2    !< EI equilibrium-kinetics analytical mixing option 2 (ideal)
        procedure :: reactive_mixing_iter_EI_eq_kin_anal_opt2    !< EI equilibrium-kinetics analytical option 2
        procedure :: reactive_mixing_iter_RKF45_eq_kin_anal_opt2 !< RKF45 equilibrium-kinetics analytical option 2 (same signature as the EI driver)
        procedure :: reactive_mixing_iter_EI_eq_anal_ideal  !< Alternative EfI option (commented)
        procedure :: reactive_mixing_iter_EI_kin_anal_ideal_opt2       !< EI kinetics analytical option 2
        procedure :: compute_u_mix                                !< Compute tilde molalities for transport
        procedure :: compute_react_term_EE_eq_kin                   !< Explicit-Euler equilibrium-kinetics reaction iteration
        procedure :: compute_react_term_EE_kin                      !< Explicit-Euler kinetics reaction iteration
        procedure :: compute_dfk_dc_EI_ideal                          !< Compute kinetics Jacobian (EfI method, arithmetic form; not used by the log-space solvers)
        procedure :: compute_dfk_dp_EI_ideal                    !< Log-space Newton Jacobian dfk/dp, p = log10 c1 (EfI, ideal)
        procedure :: compute_dfk_dp_EI                          !< Log-space Newton Jacobian dfk/dp, p = log10 c1 (EfI, non-ideal)
        procedure :: Newton_EI_eq_kin_anal_ideal_opt2         !< Newton-Raphson EI equilibrium-kinetics option 2
        procedure :: Newton_EI_eq_kin_anal_opt2                     !< Solve Newton-Raphson EI equilibrium-kinetics option 2
        procedure :: RKF45_eq_kin_anal_opt2                         !< Adaptive Runge-Kutta-Fehlberg 4(5) integration of the kinetic source, equilibrium-kinetics option 2
        procedure :: Newton_EI_kin_anal_ideal_opt2            !< Newton-Raphson EI kinetics option 2
        procedure :: Newton_EI_eq_ideal                     !< Newton-Raphson EI equilibrium-ideal
        procedure :: mix_within_zone_core                           !< Core procedure for mixing within a zone
        procedure :: mix_within_zone_new_rz
    !> Update procedures for concentration and state variables
        procedure :: update_conc_aq_prim_species                    !< Update aqueous primary species concentrations
        procedure :: update_conc_prim_species                       !< Update all primary species concentrations
        procedure :: update_conc_aq_species                         !< Update all aqueous species concentrations
        procedure :: update_state_variables               !< Update old attributes of the aqueous chemistry object
        procedure :: change_spec_alg_aq_chem                        !< Change speciation algorithm for aqueous chemistry
        procedure :: update_conc_nc                                 !< Update non-constant concentrations
        procedure :: update_attributes_dep_solid_chem                     !< Update attributes dependent on solid chemistry
        procedure :: update_attributes                        !< Update all attributes of the aqueous chemistry object
    !> Check procedures for validation and quality control
        procedure :: check_conc_aq_var_act_species                  !< Check aqueous variable activity species concentrations
        procedure :: check_conc_var_act_species                     !< Check all variable activity species concentrations
        procedure :: check_act_aq_species                           !< Check aqueous species activities for consistency
        procedure :: check_dc2v_dc1_aq                             !< Check aqueous Jacobian matrix elements
        procedure :: check_dc2v_dc1                                !< Check non-constant Jacobian matrix
        procedure :: check_dc2_dc1                                  !< Check general Jacobian matrix
    !> Initialise procedures for setting up aqueous chemistry states
        procedure :: initialise_conc_incr_coeff                     !< Initialize concentrations using incremental coefficients
        procedure :: initialise_conc_anal                           !< Initialize concentrations analytically
        procedure :: initialise_conc_anal_ideal                     !< Initialize concentrations analytically (ideal)
        procedure :: initialise_conc_anal_exch                      !< Initialize concentrations with ion exchange
        procedure :: initialise_c1_aq_iterative_method              !< Initialize primary aqueous concentrations iteratively
    !> Rearrange procedures for data structure manipulation
        procedure :: rearrange_state_vars                           !< Rearrange state variables for solver compatibility
    !> Assign procedures for object assignment and copying
        procedure :: copy_aq_chem                       !< Copy aqueous chemistry object to another
        procedure :: copy_aq_chem_Lagr                       !< Copy aqueous chemistry object to another
        procedure :: copy_aqueous_attributes                        !< Copy only the aqueous attributes of the aqueous chemistry object
        procedure :: modify_mix_ratios_reacts
        procedure :: equilibrate_with_phases                        !< Equilibrate this water with its reactive zone's minerals/gases (PHREEQC EQUILIBRIUM_PHASES style); returns per-primary-species mole transfer
    end type                                                                !< End of aqueous_chemistry_c type definition

    !> [ZONE-DBG] per-cell zone-mixing debug trace (compute_c_mix_by_zones / mix_within_zone).
    !> A solver opens within_zone_dbg_unit once (status='replace') before its time loop, sets
    !> within_zone_dbg_step each step and within_zone_dbg_cells (the target-water GLOBAL indices to
    !> trace), and closes the unit afterwards; the two subroutines only write to it (no per-call
    !> open/close). within_zone_dbg is the master enable; a block is emitted only for the selected
    !> cells during the first within_zone_dbg_max_step steps. within_zone_dbg_active is armed per
    !> selected cell by compute_c_mix_by_zones so mix_within_zone knows when to write its per-zone
    !> rows. within_zone_dbg_unit=-1 means "not open".
    public :: within_zone_dbg, within_zone_dbg_unit, within_zone_dbg_step, within_zone_dbg_cells
    public :: sim_num_time_steps, spec_fail_print_max_steps
    logical :: within_zone_dbg = .true.                   !< master enable for the zone-mixing trace
    integer :: within_zone_dbg_unit = -1                 !< output unit (opened by the solver; -1=not open)
    integer :: within_zone_dbg_step = 0                  !< current time step (set by the solver each step)
    integer :: sim_num_time_steps = 0                    !< total number of time steps of the simulation (set by the solver before the time loop; 0 = unknown)
    integer, parameter :: spec_fail_print_max_steps = 100 !< the speciation iteration-budget failure is printed only when sim_num_time_steps <= this value
    integer, parameter :: within_zone_dbg_max_step = 10  !< only trace the first this-many steps
    integer, allocatable :: within_zone_dbg_cells(:)     !< target-water GLOBAL indices to trace
    logical :: within_zone_dbg_active = .false.          !< armed by compute_c_mix_by_zones per selected cell

    !> Within-zone sub-mixture of the cross-zone mixing, written to the <root>_zone_waters.out file.
    !> chemistry_c%zone_waters(i,z) holds one per domain target water i and reactive zone z (its scratch
    !> solid chemistry is chemistry_c%zone_solids(i,z), linked once by associate_zone_waters_zone_solids).
    !> cross_zone_mix_core fills the zone waters of the target it mixes: it first resets lambda_c on all of
    !> them and then builds (lambda_c>0) only the zones of the target's mixing stencil, so at the output
    !> steps write_zone_waters (chemistry_m) selects the CURRENT step's sub-mixtures by lambda_c>0.
    type, public :: zone_water_c
        integer(kind=4) :: tar_wat_idx = 0   !< global index of the target water being mixed
        integer(kind=4) :: zone_id = 0       !< reactive zone id of the sub-mixture
        real(kind=8) :: lambda_c = 0d0       !< zone total concentration mixing ratio (its cross-zone weight); 0 = not built in the current step
        type(aqueous_chemistry_c) :: aq_chem   !< speciated sub-mixture; conc_mix holds its pre-speciation aqueous mix c_hat_z (conservative mix + limited previous-step kinetic term)
        !logical :: at_interface = .false.     !< its scratch solid chemistry, associated to the zone's reactive (and mineral) zone; holds the sub-mixture's mineral/exchange eq rates
    contains
        procedure :: set_tar_wat_idx
        procedure :: set_zone_id
        procedure :: set_lambda_c
        procedure :: set_aq_chem
        procedure :: set_zone_water_attributes
    end type zone_water_c

    !> Interface block for external subroutine signatures
    interface
        !> Compute non-constant concentrations from primary aqueous concentrations (ideal solution)
        subroutine compute_c2v_from_c1_aq_ideal(this,c2v)
            import aqueous_chemistry_c                                      !< Import the aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(out) :: c2v(:)                           !< Non-constant concentrations output array
        end subroutine
        
        !> Compute non-constant concentrations from primary concentrations (ideal solution)
        subroutine compute_c2v_from_c1_ideal(this,c1,log_act_coeffs,c2v)
            import aqueous_chemistry_c                                      !< Import the aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations input array
            real(kind=8), intent(in) :: log_act_coeffs(:) !> Logarithm (base 10) of activity coefficients of variable activity species (0 for aqueous species in ideal solution) [-]
            real(kind=8), intent(out) :: c2v(:)                           !< Non-constant concentrations output array
        end subroutine
    
        !> Compute secondary concentrations from primary aqueous concentrations (ideal solution)
        subroutine compute_c2_from_c1_aq_ideal(this,c2)
            import aqueous_chemistry_c                                      !< Import the aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(out) :: c2(:)                             !< Secondary concentrations output array
        end subroutine 
        
        !> Compute secondary concentrations from primary activities (ideal solution)
        subroutine compute_c2_from_c1_ideal(this,c1,log_act_coeffs,c2)
            import aqueous_chemistry_c                                      !< Import the aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary activities input array
            real(kind=8), intent(in) :: log_act_coeffs(:)                  !< Log activity coefficients for secondary species
            real(kind=8), intent(out) :: c2(:)                             !< Secondary concentrations output array
        end subroutine
    
        !> Compute non-constant concentrations iteratively using Picard method
        subroutine compute_c2v_from_c1_Picard(this,c1,c2v_ig,c2v,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import the aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations input
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for secondary variable activity concentrations
            real(kind=8), intent(out) :: c2v(:)                           !< Secondary variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
    
        !> Solve reactive mixing iteratively for transport-chemistry coupling
        subroutine solve_reactive_mixing_iter(this,c1_old,mixing_ratios_conc,conc_old,porosity,Delta_t,solver)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                              !< Primary species concentrations from previous time steps
            real(kind=8), intent(in) :: mixing_ratios_conc(:)               !< Mixing ratios array for transport-chemistry coupling
            real(kind=8), intent(in) :: conc_old(:,:)                       !< All concentrations from previous time step
            real(kind=8), intent(in) :: porosity                            !< Porosity for volume fraction calculations
            real(kind=8), intent(in) :: Delta_t                             !< Time step size for temporal discretization
            external :: solver                                              !< External linear solver subroutine
        end subroutine
        
        !> Compute non-constant concentrations from aqueous primary using Picard iteration
        subroutine compute_c2v_from_c1_aq_Picard(this,c2v_ig,c2v,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for secondary variable activity concentrations
            real(kind=8), intent(out) :: c2v(:)                           !< Secondary variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine 
        
        !> Compute secondary concentrations from aqueous primary using Picard iteration
        subroutine compute_c2_from_c1_aq_Picard(this,c2_ig,c2,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2_ig(:)                           !< Initial guess for secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: c2(:)                             !< Secondary concentrations output (dim=n_eq)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Compute secondary concentrations from primary using Picard iteration
        subroutine compute_c2_from_c1_Picard(this,c1,c2_ig,c2,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations input (dim=n_p)
            real(kind=8), intent(in) :: c2_ig(:)                           !< Initial guess for secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: c2(:)                             !< Secondary concentrations output (dim=n_eq)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Compute concentrations from molalities using Newton-Raphson method
        subroutine compute_c_nc_from_u_Newton(this,c1_ig,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_ig(:)                               !< Initial guess for primary concentrations
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for secondary variable activity concentrations
            real(kind=8), intent(in) :: conc_comp(:)                       !< Component concentrations constraint
            real(kind=8), intent(out) :: conc_nc(:)                        !< Variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Compute concentrations from aqueous components using Newton-Raphson method
        subroutine compute_c_nc_from_u_aq_Newton(this,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for secondary variable activity concentrations
            real(kind=8), intent(in) :: conc_comp(:)                       !< Component concentrations constraint
            real(kind=8), intent(out) :: conc_nc(:)                        !< Variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Compute concentrations from aqueous components (ideal solution, Newton-Raphson)
        subroutine compute_c_nc_from_u_aq_Newton_ideal(this,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_comp(:)                       !< Component concentrations constraint
            real(kind=8), intent(out) :: conc_nc(:)                        !< Variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Compute concentrations from components (ideal solution, Newton-Raphson)
        subroutine compute_c_nc_from_u_Newton_ideal(this,c1_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_ig(:)                           !< Initial guess for primary concentrations
            real(kind=8), intent(in) :: conc_comp(:)                       !< Component concentrations constraint
            real(kind=8), intent(out) :: conc_nc(:)                        !< Variable activity concentrations (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
                
        !> Initialize concentrations using incremental coefficients method
        subroutine initialise_conc_incr_coeff(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            integer(kind=4), intent(in) :: icon(:)                         !< Icon constraint indices array
            integer(kind=4), intent(in) :: n_icon(:)                       !< Number of icon constraints per species
            integer(kind=4), intent(in) :: indices_constrains(:,:)             !< Constraint indices for each species
            real(kind=8), intent(in) :: ctot(:)                            !< Total component concentrations constraint
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine 
        
        !> Initialize concentrations analytically with constraint handling
        subroutine initialise_conc_anal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            integer(kind=4), intent(in) :: icon(:)                         !< Icon constraint indices array
            integer(kind=4), intent(in) :: n_icon(:)                       !< Number of icon constraints per species
            integer(kind=4), intent(in) :: indices_constrains(:,:)         !< Constraint indices for each species
            real(kind=8), intent(in) :: ctot(:)                            !< Total component concentrations constraint
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
            !character(len=*), intent(in), optional :: path_py              !< Optional path for Python output
        end subroutine
        
        !> Initialize concentrations analytically (ideal solution assumption)
        subroutine initialise_conc_anal_ideal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            integer(kind=4), intent(in) :: icon(:)                         !< Icon constraint indices array
            integer(kind=4), intent(in) :: n_icon(:)                       !< Number of icon constraints per species
            integer(kind=4), intent(in) :: indices_constrains(:,:)         !< Constraint indices for each species
            real(kind=8), intent(in) :: ctot(:)                            !< Total component concentrations constraint
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
            !character(len=*), intent(in), optional :: path_py              !< Optional path for Python output
        end subroutine
        
        !> Initialize concentrations analytically with ion exchange reactions
        subroutine initialise_conc_anal_exch(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            import solid_chemistry_c                                        !< Import solid chemistry class for exchange
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            integer(kind=4), intent(in) :: icon(:)                         !< Icon constraint indices array
            integer(kind=4), intent(in) :: n_icon(:)                       !< Number of icon constraints per species
            integer(kind=4), intent(in) :: indices_constrains(:,:)         !< Constraint indices for each species
            real(kind=8), intent(in) :: ctot(:)                            !< Total component concentrations constraint
            !class(solid_chemistry_c), intent(inout) :: surf_chem          !< Surface chemistry object (commented)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine 
        
        
        !> Compute Jacobian matrix with constant activity coefficients
        subroutine compute_dc2v_dc1_gamma_cst(this,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Jacobian matrix with constant activity coefficients
        end subroutine
        
        !> Compute Jacobian matrix for ideal solution assumption
        subroutine compute_dc2_dc1_ideal(this,c1,c2,dc2_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations (dim=n_p)
            real(kind=8), intent(in) :: c2(:)                              !< Secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: dc2_dc1(:,:)                      !< Jacobian matrix for ideal solution
        end subroutine
        
        !> Compute non-constant Jacobian matrix for ideal solution
        subroutine compute_dc2v_dc1_ideal(this,c1,c2v,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations (dim=n_p)
            real(kind=8), intent(in) :: c2v(:)                            !< Non-constant secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Non-constant Jacobian matrix
        end subroutine
        
        !> Compute aqueous non-constant Jacobian matrix for ideal solution
        subroutine compute_dc2v_dc1_aq_ideal(this,c2v,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2v(:)                            !< Non-constant secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Non-constant Jacobian matrix for aqueous ideal solution
        end subroutine
        
        !> Compute non-constant Jacobian matrix with activity coefficient corrections
        subroutine compute_dc2v_dc1(this,c1,c2v,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations for activity calculations
            real(kind=8), intent(in) :: c2v(:)                            !< Non-constant secondary concentrations
            !real(kind=8), intent(in) :: out_prod(:,:)                      !< Outer product between d_log_gamma_nc_d_I and z_nc^2 (activity coefficient derivatives)
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Non-constant Jacobian matrix with activity corrections
        end subroutine

        !> Log-log Jacobian d(log10 c2v)/d(log10 c1) with activity coefficient corrections (the quantity
        !> compute_dc2v_dc1 computes before converting it to dc2v/dc1 = c2v L / c1); used by the log-space
        !> assembly of the Newton Jacobian (compute_dfk_dp_EI) so that no 1/c factor is ever formed
        subroutine compute_dlogc2v_dlogc1(this,c1,c2v,dlogc2v_dlogc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic), speciated at (c1,c2v)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations, speciation order (n_p)
            real(kind=8), intent(in) :: c2v(:)                            !< Secondary variable activity concentrations, speciation order (n_e)
            real(kind=8), intent(out) :: dlogc2v_dlogc1(:,:)              !< d(log10 c2v)/d(log10 c1) (n_e x n_p, already allocated)
        end subroutine
        
        !> Compute aqueous non-constant Jacobian matrix with activity coefficient corrections
        subroutine compute_dc2v_dc1_aq(this,c2v,out_prod,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2v(:)                            !< Non-constant secondary concentrations
            real(kind=8), intent(in) :: out_prod(:,:)                      !< Outer product between d_log_gamma_nc_d_I and z_nc^2 (activity derivatives)
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Aqueous non-constant Jacobian with activity corrections
        end subroutine
                
        !> Compute general Jacobian matrix for secondary species with activity corrections
        subroutine compute_dc2_dc1(this,out_prod,c1,c2,dc2_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: out_prod(:,:)                      !< Outer product between d_log_gamma_d_I and z^2 (activity coefficient derivatives)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations (dim=n_p)
            real(kind=8), intent(in) :: c2(:)                              !< Secondary concentrations (dim=n_eq)
            real(kind=8), intent(out) :: dc2_dc1(:,:)                      !< General Jacobian matrix with activity corrections
        end subroutine
        
        
        !> Compute non-constant concentrations explicitly for homogeneous systems
        subroutine compute_c2v_from_c1_expl_homog(this)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic) - explicit computation for homogeneous phase
        end subroutine 
        
        !> Compute non-constant concentrations using Picard iteration for homogeneous systems
        subroutine compute_c2v_from_c1_Picard_homog(this,tolerance,niter_max,niter)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: tolerance                              !< Convergence tolerance for Picard iteration
            integer(kind=4), intent(in) :: niter_max                       !< Maximum number of iterations allowed
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
        end subroutine 
        
        !> Compute concentrations from molalities using Newton-Raphson for homogeneous systems
        subroutine compute_c_nc_from_u_Newton_homog(this,tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: tolerance                          !< Absolute convergence tolerance
            real(kind=8), intent(in) :: rel_tolerance                      !< Relative convergence tolerance
            real(kind=8), intent(in) :: control_factor                     !< Newton step control factor for stability
            integer(kind=4), intent(in) :: niter_max                       !< Maximum number of iterations allowed
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine 
        
        !> Compute residual vector for homogeneous equilibrium chemistry
        subroutine compute_res_spec_homog(this,residual)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (read-only)
            real(kind=8), intent(out) :: residual(:)                       !< Residual vector for equilibrium equations
        end subroutine
        
        !> Compute Jacobian matrix with constant activity coefficients for homogeneous systems
        subroutine compute_dc2v_dc1_gamma_cst_homog(this,tol,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (read-only)
            real(kind=8), intent(in) :: tol                                !< Tolerance for numerical differentiation
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Jacobian matrix with constant activity coefficients
        end subroutine
        
        !> Compute full Jacobian matrix for homogeneous systems with activity corrections
        subroutine compute_dc2v_dc1_homog(this,tol,dc2v_dc1)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (read-only)
            real(kind=8), intent(in) :: tol                                !< Tolerance for numerical differentiation
            real(kind=8), intent(out) :: dc2v_dc1(:,:)                    !< Full Jacobian matrix with activity coefficient derivatives
        end subroutine 
        
        subroutine reactive_mixing_iter_EI_eq_anal_ideal(this,c1_old,c_hat,mix_ratio_r,&
            Delta_t,theta,conc_nc,conc_comp,c1_downstream)    !< Update aqueous species concentrations (commented)
            import aqueous_chemistry_c                                   !< Import aqueous chemistry class
            implicit none                                                !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                         !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:) !> aqueous primary concentrations at previous time step
            real(kind=8), intent(in) :: c_hat(:) !> variable activity species concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r !> reaction mixing ratio of this target water
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations (already allocated)
            real(kind=8), intent(inout) :: conc_comp(:) !> at this target
            real(kind=8), intent(in), optional :: c1_downstream(:) !> primary concentrations of closest downstream target water
        end subroutine
        
    ! !> Updates concentration aqueous primary species in iterative method (commented)
        subroutine Newton_EI_eq_ideal(this,u_hat,conc_nc,niter,CV_flag)       !< Update aqueous primary species concentrations
            import aqueous_chemistry_c                               !< Import aqueous chemistry class
            implicit none                                            !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                      !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: u_hat(:) !> Component concentrations after mixing [C] (size = n_comp)
            real(kind=8), intent(inout) :: conc_nc(:) !> Variable activity species concentrations c [C] (size = n_v, pre-allocated, updated in Newton loop)
            integer(kind=4), intent(out) :: niter !> Number of Newton iterations performed [-]
            logical, intent(out) :: CV_flag !> Convergence flag: FALSE = no convergence, TRUE = converged
        end subroutine
        
        
        !> Compute logarithmic activity coefficients for all aqueous species
        subroutine compute_log_act_coeffs(this)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic) - computes activity coefficients
        end subroutine
        
        !> Compute Jacobian matrix of logarithmic activity coefficients with respect to concentrations
        subroutine compute_log_Jacobian_act_coeffs(this,d_log_gamma_d_log_c)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(out) :: d_log_gamma_d_log_c(:,:)          !< Jacobian: d(log gamma)/d(log c) for activity coefficient sensitivity
        end subroutine
        
        !> Component mixing iteration with reactive transport coupling
        subroutine mixing_iter_comp(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for non-constant secondary concentrations
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            !real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                       !< Mixing ratio for kinetic reactions
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
        end subroutine

        !> Component mixing iteration with lumped kinetics approach
        subroutine mixing_iter_comp_lump(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                    !< Mixing ratio for kinetic reactions
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
        end subroutine
        
        !> Component mixing iteration with ideal solution assumption
        !> Same signature as mixing_iter_comp so both can share a procedure(mixing_iter_comp) pointer
        subroutine mixing_iter_comp_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            ! real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                       !< Mixing ratio for kinetic reactions
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
            !real(kind=8), intent(in), optional :: Delta_t                 !< Optional time step (commented)
        end subroutine

        !> Speciate a water with its reactive zone (ideal): copy of mixing_iter_comp_ideal WITHOUT the
        !> unchanged-primaries early-exit (Newton always runs). Used by mix_within_zone on displaced waters.
        subroutine equilibrate_water_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c_hat(:)
            real(kind=8), intent(in) :: mix_ratio_r
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(inout) :: conc_nc(:)
            real(kind=8), intent(inout) :: conc_comp(:)
        end subroutine

        !> Speciate a water with its reactive zone (non-ideal): copy of mixing_iter_comp WITHOUT the
        !> unchanged-primaries early-exit (Newton always runs). Used by mix_within_zone on displaced waters.
        subroutine equilibrate_water(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c_hat(:)
            real(kind=8), intent(in) :: mix_ratio_r
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(inout) :: conc_nc(:)
            real(kind=8), intent(inout) :: conc_comp(:)
        end subroutine

        !> Component mixing iteration with ideal solution and lumped kinetics
        subroutine mixing_iter_comp_ideal_lump(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal lumped)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            !real(kind=8), intent(in) :: rk_tilde(:)                       !< Kinetic reaction rates (commented - lumped approach)
            real(kind=8), intent(in) :: mix_ratio_r                      !< Mixing ratio (commented - lumped approach)
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
            !real(kind=8), intent(in), optional :: Delta_t                 !< Optional time step (commented)
        end subroutine
        
        !> Component mixing iteration with ion exchange reactions
        subroutine mixing_iter_comp_exch(this,c1_old,c2v_ig,c_mix,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for non-constant secondary concentrations
            real(kind=8), intent(in) :: c_mix(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_Rk                       !< Mixing ratio for kinetic reactions
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            !real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
        end subroutine

        !> Component mixing iteration with ion exchange and lumped kinetics
        subroutine mixing_iter_comp_exch_lump(this,c1_old,c2v_ig,c_mix,Delta_t,conc_nc)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for non-constant secondary concentrations
            real(kind=8), intent(in) :: c_mix(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            !real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
        end subroutine
        
        !> Component mixing iteration with ion exchange and ideal solution assumption
        subroutine mixing_iter_comp_exch_ideal(this,c1_old,c_mix,conc_nc,porosity,Delta_t)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_mix(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            !real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
            real(kind=8), intent(in), optional :: porosity                 !< Optional porosity for volume fraction calculations
            real(kind=8), intent(in), optional :: Delta_t                  !< Optional time step size for temporal discretization
        end subroutine

        !> Component mixing iteration with ion exchange, ideal solution, and lumped kinetics
        subroutine mixing_iter_comp_exch_ideal_lump(this,c1_old,c_mix,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal lumped)
            real(kind=8), intent(in) :: c_mix(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            !real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
        end subroutine
        
        !> Species-level mixing iteration for detailed chemical speciationm
        subroutine mixing_iter_species(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations before previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            ! real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                    !< Mixing ratio for reactions (modified)
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
        end subroutine
            
        
        !> Mixing iteration using Explicit-Implicit (EI) method for aqueous chemistry
        subroutine mixing_iter_EI_aq_chem(this,c_mix,rk_mat,porosity,Delta_t)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c_mix(:,:)                       !< Transport-advected concentration matrix (spatial nodes)
            real(kind=8), intent(in), optional :: rk_mat(:,:)              !< Optional kinetic reaction rate matrix
            real(kind=8), intent(in), optional :: porosity                 !< Optional porosity storage matrix
            real(kind=8), intent(in), optional :: Delta_t                  !< Optional time step size
        end subroutine
        
        !> Mixing iteration using Explicit fully Implicit (EfI) method for aqueous chemistry
        subroutine mixing_iter_EfI_aq_chem(this,porosity,Delta_t,tolerance,rel_tolerance,control_factor,niter_max)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            !class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)    !< All waters have same aqueous species and components (commented)
            !real(kind=8), intent(in) :: mixing_ratios_conc(:)             !< Mixing ratios A^(-T)_j (commented)
            real(kind=8), intent(in), optional :: porosity                 !< Optional porosity storage matrix
            real(kind=8), intent(in), optional :: Delta_t                  !< Optional time step size
            real(kind=8), intent(in), optional :: tolerance                !< Optional absolute tolerance for Newton convergence
            real(kind=8), intent(in), optional :: rel_tolerance            !< Optional relative tolerance for Newton/Picard convergence
            real(kind=8), intent(in), optional :: control_factor           !< Optional control factor for Delta_c1_j in Newton method
            integer(kind=4), intent(in), optional :: niter_max             !< Optional maximum number of iterations
        end subroutine
        
        !> Water mixing iteration using Explicit-Explicit (EE) method with equilibrium and kinetics
        subroutine reactive_mixing_iter_EE_eq_kin(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for non-constant secondary concentrations
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            !real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                       !< Mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                        !< Non-constant concentrations output
        end subroutine
        
        !> Water mixing iteration using EE method with equilibrium and kinetics (ideal solution)
        subroutine reactive_mixing_iter_EE_eq_kin_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            !real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                       !< Mixing ratio of reaction amounts in this target
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity storage matrix (commented)
            
        end subroutine
                
        !> Water mixing iteration using EE method with kinetics only (ideal solution)
        subroutine reactive_mixing_iter_EE_kin(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            !real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                       !< Mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity storage matrix (commented)
            
        end subroutine
        
        !> Water mixing iteration using EfI method with equilibrium and kinetics (analytical, lumped)
        subroutine reactive_mixing_iter_EfI_eq_kin_anal_lump(this,c1_old,c2v_ig,c_mix,Delta_t,conc_nc)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            !import time_discr_c                                            !< Time discretization module (commented)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            !class(aqueous_chemistry_c), intent(in) :: this_old            !< Aqueous chemistry object at previous time step (commented)
            !integer(kind=4), intent(in) :: ind                            !< Index of object "this" in mixing ratios array (commented)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c2v_ig(:)                         !< Initial guess for non-constant secondary concentrations
            real(kind=8), intent(in) :: c_mix(:)                         !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
        end subroutine
        
        !> Water mixing iteration using EI method with equilibrium and kinetics (analytical, ideal, option 2)
        subroutine reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2(this,c1_old,c_hat,&
            mix_ratio_r,Delta_t,theta,conc_nc,conc_comp,c1_downstream)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                            !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c_hat(:)                          !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                        !< Reaction mixing ratio of this target water
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            real(kind=8), intent(in), optional :: c1_downstream(:)          !< Primary concentrations of closest downstream target water
        end subroutine
        
        !> Water mixing iteration using EI method with equilibrium and kinetics (analytical, non-ideal, option 2)
        subroutine reactive_mixing_iter_EI_eq_kin_anal_opt2(this,c1_old,c_hat,&
            mix_ratio_r,Delta_t,theta,conc_nc,conc_comp,c1_downstream)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                            !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c_hat(:)                          !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                        !< Reaction mixing ratio of this target water
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            real(kind=8), intent(in), optional :: c1_downstream(:)          !< Primary concentrations of closest downstream target water
        end subroutine

        !> Water mixing iteration with adaptive RKF45 kinetics, equilibrium and kinetics (analytical, non-ideal, option 2)
        subroutine reactive_mixing_iter_RKF45_eq_kin_anal_opt2(this,c1_old,c_hat,&
            mix_ratio_r,Delta_t,theta,conc_nc,conc_comp,c1_downstream)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                            !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c_hat(:)                          !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                        !< Reaction mixing ratio of this target water
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components
            real(kind=8), intent(in), optional :: c1_downstream(:)          !< Primary concentrations of closest downstream target water
        end subroutine

        subroutine reactive_mixing_iter_EI_kin_anal_ideal_opt2(this,c1_old,c_hat,&
            mix_ratio_r,Delta_t,&
            theta,conc_nc,conc_comp,c1_downstream)
            import aqueous_chemistry_c
            !import diag_matrix_c
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2v_ig(:)
            real(kind=8), intent(in) :: c_hat(:)
            !real(kind=8), intent(in) :: rk_tilde_up(:)
            !real(kind=8), intent(in) :: rk_tilde_down(:)
            real(kind=8), intent(in) :: mix_ratio_r !> reaction mixing ratio of this target water
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:)
            real(kind=8), intent(inout) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            real(kind=8), intent(in), optional :: c1_downstream(:) !> primary concentrations of closest downstream target water
        end subroutine

        subroutine reactive_mixing_iter_EI_kin_anal_ideal_opt2_bis(this,c1_old,&
            c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
           import aqueous_chemistry_c
           import diag_matrix_c 
           !import time_discr_c
           implicit none
           class(aqueous_chemistry_c) :: this
           real(kind=8), intent(in) :: c1_old(:)
           !real(kind=8), intent(in) :: c2v_ig(:)
           real(kind=8), intent(in) :: c_hat(:)
           real(kind=8), intent(in) :: mix_ratio_r !> reaction mixing ratio of this target water
           real(kind=8), intent(in) :: Delta_t !> time step
           real(kind=8), intent(in) :: theta !> time weighting factor
           real(kind=8), intent(inout) :: conc_nc(:)
           real(kind=8), intent(inout) :: conc_comp(:)
           !real(kind=8), intent(in), optional :: porosity
        end subroutine

        subroutine compute_react_term_EE_eq_kin(this,Delta_t,lambda_r,conc_comp_react)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: lambda_r !> kinetic reaction rate contributions after mixing
            real(kind=8), intent(out) :: conc_comp_react(:) !> must be allocated
        end subroutine
        
        subroutine compute_react_term_EE_eq_kin_lump(this,Delta_t,conc_comp_react)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
            real(kind=8), intent(out) :: conc_comp_react(:) !> must be allocated
        end subroutine

        subroutine compute_react_term_EE_kin(this,Delta_t,lambda_r,conc_react)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: lambda_r !> kinetic reaction rate contributions after mixing
            real(kind=8), intent(out) :: conc_react(:) !> reaction part of variable activity concentrations (must be already allocated)
        end subroutine

        
        !> Jacobian of the Newton residual of the reactive mixing with respect to p = log10 c1 (log-space
        !> assembly from the log-form kinetic Jacobian drk_dlogc of compute_rk_Jac_rk_anal; ideal solution)
        subroutine compute_dfk_dp_EI_ideal(this,c1,c2v,drk_dlogc,Delta_t,theta,mix_ratio_r,dfk_dp)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c1(:)
            real(kind=8), intent(in) :: c2v(:)
            real(kind=8), intent(in) :: drk_dlogc(:,:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: mix_ratio_r
            real(kind=8), intent(out) :: dfk_dp(:,:)
        end subroutine

        !> Same as compute_dfk_dp_EI_ideal for a NON-ideal solution (activity coefficients depending on the
        !> ionic strength: the secondary-primary chain rule is compute_dlogc2v_dlogc1)
        subroutine compute_dfk_dp_EI(this,c1,c2v,drk_dlogc,Delta_t,theta,mix_ratio_r,dfk_dp)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c1(:)
            real(kind=8), intent(in) :: c2v(:)
            real(kind=8), intent(in) :: drk_dlogc(:,:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: mix_ratio_r
            real(kind=8), intent(out) :: dfk_dp(:,:)
        end subroutine

         subroutine compute_dfk_dc_EI_ideal(this,drk_dc,Delta_t,theta,mix_ratio_r,dfk_dc)
             import aqueous_chemistry_c
             implicit none
             class(aqueous_chemistry_c), intent(in) :: this
             real(kind=8), intent(in) :: drk_dc(:,:)
             real(kind=8), intent(in) :: Delta_t
             real(kind=8), intent(in) :: theta
             real(kind=8), intent(in) :: mix_ratio_r
             real(kind=8), intent(out) :: dfk_dc(:,:)
         end subroutine
        
      
        subroutine Newton_EI_eq_kin_anal_ideal_opt2(this,u_hat,mix_ratio_r_new,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine

        subroutine Newton_EI_eq_kin_anal_opt2(this,u_hat,mix_ratio_r_new,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine

        subroutine RKF45_eq_kin_anal_opt2(this,u_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: u_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of accepted RKF45 sub-steps
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        
        subroutine Newton_EI_kin_anal_ideal_opt2(this,c_hat,&
            mix_ratio_r,Delta_t,theta,conc_nc,niter,&
            CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
            end subroutine

            subroutine Newton_EI_kin_anal_ideal_opt2_bis(this,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,&
                CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: c_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio of kinetic reaction rate in this target
            !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
            end subroutine
            
      
        subroutine compute_re_mean_tilde_aq_chem(this,this_old,mixing_ratios,mixing_waters,Delta_t,re_mean_tilde)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            class(aqueous_chemistry_c), intent(in) :: this_old              !> nombre muy cutre
            real(kind=8), intent(in) :: mixing_ratios(:)
            class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(out) :: re_mean_tilde(:)
        end subroutine
        
        
        subroutine compute_Re(this,cv_hat,Delta_t,lambda_r,theta_r)
            import aqueous_chemistry_c
            implicit none  
            class(aqueous_chemistry_c) :: this !> object
            real(kind=8), intent(in) :: cv_hat(:) !> concentrations of all variable activity species after mixing
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: lambda_r !> mixing ratio factor for reaction rates
            real(kind=8), intent(in) :: theta_r !> reaction time weighting factor (weights only the kinetic term)
        end subroutine
        
        subroutine compute_Re_lump(this,c2v_mix,Delta_t,theta)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2v_mix(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: porosity
        end subroutine

        subroutine compute_Re_rk_lump(this,c2v_mix,Delta_t,theta)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2v_mix(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: porosity
        end subroutine
        
        subroutine compute_rk_eq_heterog(this)!,rk)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
        end subroutine
        
        subroutine compute_rxn_rate_kin(this,rk_new)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: rk_new(:)
        end subroutine
        
        subroutine compute_Jacobian_rk_anal(this,drk_dlogc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dlogc(:,:) !> Jacobian of the kinetic rates in LOG form, d r_k / d log10 c (chem_syst column order)
        end subroutine
        
        subroutine compute_rk_Jac_rk_anal(this,rk_new,drk_dlogc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: rk_new(:) !>  (must be allocated)
            real(kind=8), intent(out) :: drk_dlogc(:,:) !> Jacobian of the kinetic rates in LOG form, d r_k / d log10 c = ln10 c d r_k / d c (chem_syst column order)
        end subroutine
        
        subroutine compute_rk_Jac_rk_incr_coeff(this,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:,:)
        end subroutine
        
   
        function compute_d_log_gamma_d_I_aq_chem(this) result(d_log_gamma_d_I)
            import aqueous_chemistry_c
            import params_aq_sol_s
            implicit none
            class(aqueous_chemistry_c) :: this !> object
            real(kind=8), allocatable :: d_log_gamma_d_I(:) !> Derivative of logÃ¢â€šÂÃ¢â€šâ‚¬(activity coefficients) w.r.t. ionic strength: Ã¢Ë†â€š(log ÃŽÂ³)/Ã¢Ë†â€šI (must be already allocated) [MÃ¢ÂÂ»Ã‚Â¹]
        end function
        
        subroutine compute_res_Jac_res_anal_exch(this,conc,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,&
            log_Jacobian_act_coeffs,CEC,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            real(kind=8), intent(in) :: CEC
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_anal(this,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,log_Jacobian_act_coeffs,res,&
            Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_anal_ideal(this,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_incr_coeff(this,c2_init,indices_icon,n_icon,indices_constrains,ctot,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2_init(:)
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_init(this,indices_icon,n_icon,indices_constrains,ctot,res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            !> Pre-process
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
        end subroutine
        
        subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,num_gas_zones,model,Jac_opt,unit,niter,CV_flag,gas_species_chem)
            import aqueous_chemistry_c
            !import chem_system_c
            import solid_chemistry_c
            import gas_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
            integer(kind=4), intent(in) :: num_cstr !> number of constrains
            integer(kind=4), intent(in) :: num_gas_zones !> number of gas zones
            integer(kind=4), intent(in) :: model !> activity coefficients model
            integer(kind=4), intent(in) :: Jac_opt !> 0: incremental coeffficinets, 1: analtical
            integer(kind=4), intent(in) :: unit !> file unit
            integer(kind=4), intent(out) :: niter !> number of iterations in Newton method
            logical, intent(out) :: CV_flag !> TRUE if Newton method converges, FALSE otherwise
            !type(solid_chemistry_c), intent(in), optional :: surf_chem
            type(gas_chemistry_c), intent(in), optional :: gas_species_chem
        end subroutine
        
          
    end interface                                                           !< End of final interface block

    !> Guarded copy of an allocatable array component: dst takes the value of src when src is
    !> allocated and is deallocated otherwise. A plain "dst = src" with an UNALLOCATED src is not a
    !> valid reference (gfortran reads the undefined bounds of the empty descriptor and (re)allocates
    !> dst to a garbage size): the deallocate-on-unallocated-source rule of intrinsic assignment only
    !> holds for allocatable COMPONENTS copied through a derived-type assignment.
    interface copy_alloc
        module procedure copy_alloc_r1, copy_alloc_r2, copy_alloc_i1
    end interface copy_alloc
    private :: copy_alloc_r1, copy_alloc_r2, copy_alloc_i1
    
    contains                                                                !< Beginning of implementation section

        subroutine copy_alloc_r1(dst,src)
            real(kind=8), allocatable, intent(inout) :: dst(:)
            real(kind=8), allocatable, intent(in) :: src(:)
            if (allocated(src)) then
                dst=src
            else if (allocated(dst)) then
                deallocate(dst)
            end if
        end subroutine copy_alloc_r1

        subroutine copy_alloc_r2(dst,src)
            real(kind=8), allocatable, intent(inout) :: dst(:,:)
            real(kind=8), allocatable, intent(in) :: src(:,:)
            if (allocated(src)) then
                dst=src
            else if (allocated(dst)) then
                deallocate(dst)
            end if
        end subroutine copy_alloc_r2

        subroutine copy_alloc_i1(dst,src)
            integer(kind=4), allocatable, intent(inout) :: dst(:)
            integer(kind=4), allocatable, intent(in) :: src(:)
            if (allocated(src)) then
                dst=src
            else if (allocated(dst)) then
                deallocate(dst)
            end if
        end subroutine copy_alloc_i1
    
     
        !> Allocate memory for concentrations of aqueous species after conservative mixing of waters
        subroutine allocate_conc_mix(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%conc_mix)) deallocate(this%conc_mix) !< Deallocate if already allocated
            allocate(this%conc_mix(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))         !< Allocate current concentration array for all aqueous species
        end subroutine

        !> Allocate memory for the aqueous species concentrations of the mixed waters entering the Newton speciation
        subroutine allocate_conc_hat(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%conc_hat)) deallocate(this%conc_hat) !< Deallocate if already allocated
            allocate(this%conc_hat(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !< [aq prim | aq sec var act] of the water's reactive zone
            this%conc_hat=0d0
        end subroutine

        !> Store the aqueous blocks of the Newton speciation input c_hat (full variable-activity layout of
        !> the water's reactive zone: [all primary | aq secondary var act | surface complexes | gases])
        !> in conc_hat ([aq prim | aq sec var act]). Allocates conc_hat if needed.
        subroutine set_conc_hat(this,c_hat)
            implicit none
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c_hat(:)                           !< speciation input, size >= num_prim_species+num_aq_sec_var_act_species
            integer(kind=4) :: n_p,n_aq_p,n_aq_s
            associate(sa => this%solid_chemistry%reactive_zone%speciation_alg)
            n_p=sa%num_prim_species
            n_aq_p=sa%num_aq_prim_species
            n_aq_s=sa%num_aq_sec_var_act_species
            end associate
            if (size(c_hat)<n_p+n_aq_s) error stop "set_conc_hat: c_hat shorter than the primary + aqueous secondary blocks"
            if (.not. allocated(this%conc_hat)) then
                call this%allocate_conc_hat()
            else if (size(this%conc_hat)/=n_aq_p+n_aq_s) then
                call this%allocate_conc_hat()
            end if
            this%conc_hat(1:n_aq_p)=c_hat(1:n_aq_p)
            this%conc_hat(n_aq_p+1:n_aq_p+n_aq_s)=c_hat(n_p+1:n_p+n_aq_s)
        end subroutine

        !> Allocate memory for aqueous species concentrations
        subroutine allocate_conc_aq_species(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%concentrations)) deallocate(this%concentrations)      !< Deallocate if already allocated
            if (allocated(this%conc_old)) deallocate(this%conc_old)                  !< Deallocate if already allocated
            if (allocated(this%conc_old_old)) deallocate(this%conc_old_old)         !< Deallocate if already allocated
            allocate(this%concentrations(this%aq_phase%num_species),&
                this%conc_old(this%aq_phase%num_species),&
                this%conc_old_old(this%aq_phase%num_species))              !< Allocate current, previous, and previous-previous concentration arrays
        end subroutine

        !> Allocate memory for contributions to aqueous species concentrations from kinetic and equilibrium reaction amounts
        subroutine allocate_contrib_rxn_amnt(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%contrib_rxn_amnt_kin)) deallocate(this%contrib_rxn_amnt_kin) !< Deallocate if already allocated
            if (allocated(this%contrib_rxn_amnt_eq)) deallocate(this%contrib_rxn_amnt_eq) !< Deallocate if already allocated
            allocate(this%contrib_rxn_amnt_kin(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),&
                this%contrib_rxn_amnt_eq(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))              !< Allocate the kinetic and equilibrium contribution arrays (variable activity species)
            !> Zero contributions by default: the kinetic one is only assigned by the kinetic iterations (never in an
            !> equilibrium-only system) and the equilibrium one only by compute_Re, so they are otherwise read
            !> (e.g. by the [CELL-DBG] trace) as uninitialised memory
            this%contrib_rxn_amnt_kin=0d0
            this%contrib_rxn_amnt_eq=0d0
        end subroutine
        
        !> Allocate memory for aqueous species activities
        subroutine allocate_activities_aq_species(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            allocate(this%activities(this%aq_phase%num_species))           !< Allocate activities array for all aqueous species
        end subroutine
        
        !> Compute activities for dissolved solid species using activity coefficients
        subroutine compute_activities_diss_solids(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            this%activities(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))*&
                (10**(this%log_act_coeffs(this%indices_aq_species(this%aq_phase%ind_diss_solids))))
                                                                            !< Activity = concentration Ãƒâ€” activity_coefficient (a = c Ãƒâ€” ÃŽÂ³)
        end subroutine
        
        !> Compute activities of aqueous primary species
        subroutine compute_a1(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            this%activities(this%ind_prim_species)=&
                this%concentrations(this%ind_prim_species)*&
                (10**(this%log_act_coeffs(this%ind_prim_species)))          !< Primary activities: a1 = c1 Ãƒâ€” ÃŽÂ³1
        end subroutine
        
        !> Compute activities of secondary species
        subroutine compute_a2(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            this%activities(this%ind_sec_species)=&
                this%concentrations(this%ind_sec_species)*&
                (10**(this%log_act_coeffs(this%ind_sec_species)))           !< Secondary activities: a2 = c2 Ãƒâ€” ÃŽÂ³2
        end subroutine
        
        !> Compute activities of all aqueous species (dissolved solids and water)
        !> \brief Computes activities of all aqueous species
        !!
        !! \details This subroutine updates the chemical activities of all aqueous species
        !! based on their current concentrations and activity coefficients. For aqueous
        !! solutions, activity is related to concentration via:
        !! \f[
        !! a_i = \gamma_i \cdot c_i
        !! \f]
        !! where:
        !! - \f$ a_i \f$ = activity of species i [-]
        !! - \f$ \gamma_i \f$ = activity coefficient of species i [-]
        !! - \f$ c_i \f$ = concentration of species i [C] or [M]
        !!
        !! **Algorithm**:
        !! 1. Compute activities of dissolved solid species (ions, complexes)
        !! 2. Compute water activity using Raoult's law approximation
        !!
        !! \note Activity coefficients must be computed beforehand (e.g., via Debye-HÃƒÂ¼ckel,
        !!       Davies, or Pitzer equations). This subroutine only applies them to get activities.
        !!
        !! \warning For ideal solutions (ÃŽÂ³ = 1), activity equals concentration.
        !!
        !! \see compute_activities_diss_solids, compute_act_water, compute_ionic_strength
        !!
        subroutine compute_activities_aq(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            !> \subsection compute_dissolved_activities Compute activities of dissolved species
            !! Update activities for all dissolved species (ions and complexes) in the
            !! aqueous phase using a_i = ÃŽÂ³_i * c_i
            call this%compute_activities_diss_solids()                     !< Apply activity coefficients to dissolved solid concentrations
            
            !> \subsection compute_water_activity Compute water activity
            !! Calculate water activity using Raoult's law: a_w Ã¢â€°Ë† 1 - 0.018 * ÃŽÂ£m_i
            !! where the sum is over all dissolved species molalities
            call this%compute_act_water()                                  !< Compute water activity via Raoult's law
        end subroutine
        
        !> \brief Computes activities of all species (aqueous and solid phases)
        !!
        !! \details This subroutine updates the chemical activities of all species in both
        !! the aqueous phase and the solid phase. This comprehensive update is necessary
        !! when both phases are present and coupled through reactions (e.g., mineral
        !! dissolution/precipitation). The activity calculations include:
        !! - Aqueous species: \f$ a_i = \gamma_i \cdot c_i \f$
        !! - Water: Raoult's law approximation
        !! - Solids: typically a = 1 for pure phases, or ion activity products for solid solutions
        !!
        !! **Algorithm**:
        !! 1. Compute activities of dissolved solid species (aqueous ions, complexes)
        !! 2. Compute water activity using Raoult's law
        !! 3. Compute activities of solid phase species (minerals)
        !!
        !! \note This subroutine requires solid_chemistry to be associated. Use
        !!       compute_activities_aq if only aqueous phase is present.
        !!
        !! \warning For pure minerals, activity is typically set to 1. For solid solutions
        !!          (e.g., (Ca,Mg)COÃ¢â€šÆ’), activities depend on mole fractions.
        !!
        !! \see compute_activities_aq, compute_activities_diss_solids, solid_chemistry_c::compute_activities_solids
        !!
        subroutine compute_activities(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            !> \subsection compute_aq_activities_all Compute aqueous activities
            !! Update activities for all dissolved species (ions and complexes) using
            !! the relationship a_i = ÃŽÂ³_i * c_i
            call this%compute_activities_diss_solids()                     !< Apply activity coefficients to dissolved concentrations
            
            !> \subsection compute_water_activity_all Compute water activity
            !! Calculate water activity using Raoult's law approximation:
            !! a_w Ã¢â€°Ë† 1 - 0.018 * ÃŽÂ£m_i for dilute solutions
            call this%compute_act_water()                                  !< Compute water activity via Raoult's law
            
            !> \subsection compute_solid_activities Compute solid phase activities
            !! Update activities for all mineral species. For pure phases a = 1,
            !! for solid solutions a depends on mole fractions
            call this%solid_chemistry%compute_activities_solids()          !< Compute mineral activities
        end subroutine
        
        !> Compute activity of a specific aqueous species
        subroutine compute_activity(this,ind)
        !> compute activity of an aqueous species
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            integer(kind=4), intent(in) :: ind                             !< Index of aqueous species in state variables
            this%activities(ind)=this%concentrations(ind)*(10**this%log_act_coeffs(ind))
                                                                            !< Activity calculation: a_i = c_i Ãƒâ€” ÃŽÂ³_i
        end subroutine
        
        !> Compute water activity using Raoult's law approximation
        subroutine compute_act_water(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (this%aq_phase%ind_wat>0) then                              !< If water index is defined
                this%activities(this%indices_aq_species(this%aq_phase%ind_wat))=1d0-0.018*sum(this%concentrations(&
                    this%indices_aq_species(this%aq_phase%ind_diss_solids)))
                                                                            !< First-order Raoult's law: a_w = 1 - 0.018 Ãƒâ€” ÃŽÂ£c_i
            end if
        end subroutine
                
        !> Allocate memory for logarithmic activity coefficients and their Jacobian
        subroutine allocate_log_act_coeffs_aq_chem(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            allocate(this%log_act_coeffs(this%aq_phase%num_species))       !< Allocate log activity coefficients array
            allocate(this%log_Jacobian_act_coeffs(this%aq_phase%num_species,this%aq_phase%num_species))
                                                                            !< Allocate Jacobian matrix for activity coefficient derivatives
            this%log_act_coeffs=0d0                                            !< Initialize log activity coefficients to zero (ideal solution)
            this%log_Jacobian_act_coeffs=0d0                               !< Initialize Jacobian matrix to zero
        end subroutine
        
        !> Set concentrations of all aqueous species with dimension checking
        subroutine set_conc_aq_species(this,conc)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc(:)                            !< Input concentration array
            if (size(conc)/=this%aq_phase%num_species) then                !< Check array dimension compatibility
                error stop "Dimension error in concentration of aqueous species"
            else
                this%concentrations=conc                                   !< Assign concentrations if dimensions match
            end if
        end subroutine
        
        !> Set or initialize logarithmic activity coefficients for aqueous species
        subroutine set_log_act_coeffs(this,log_act_coeffs) !> we set or initialise the logarithm of activity coefficientes of  & 
            !! aqueous species
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in), optional :: log_act_coeffs(:)        !< Optional input log activity coefficients
            integer(kind=4) :: i                                           !< Loop counter for species
            if (present(log_act_coeffs)) then                              !< If log activity coefficients are provided
                this%log_act_coeffs=log_act_coeffs                         !< Assign provided values
            else if (.not. allocated(this%log_act_coeffs)) then            !< If array not allocated
                error stop "Activity coefficients not allocated"
            else                                                           !< Default initialization
                do i=1,this%aq_phase%num_species                           !< Loop over all aqueous species
                    this%log_act_coeffs(i)=0d0                             !< Set to zero (ideal solution assumption)
                end do
            end if
        end subroutine
        
        !> Set concentrations of aqueous primary species with validation
        subroutine set_conc_aq_prim_species(this,conc_aq_prim_species)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object
            real(kind=8), intent(in) :: conc_aq_prim_species(:)            !< Aqueous primary species concentrations input
            
            integer(kind=4) :: i                                           !< Loop counter for primary species
            if (size(conc_aq_prim_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                error stop "Dimension error in concentration of primary species" !< Dimension validation error
            else                                                           !< If dimensions are correct
                do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !< Loop over primary species
                    this%concentrations(this%ind_prim_species(i))=conc_aq_prim_species(i)
                end do
            end if
        end subroutine
        
        !> Set concentrations of all primary species (aqueous and solid)
        subroutine set_conc_prim_species(this,conc_prim_species)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_prim_species(:)               !< Primary species concentrations (all phases)
            if (size(conc_prim_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) then
                error stop "Dimension error in concentration of primary species" !< Dimension validation
            end if
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1)=conc_prim_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !< Set solid primary species concentration
            end if
            this%concentrations(this%ind_prim_species)=&
                conc_prim_species(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
                                                                            !< Set aqueous primary species concentrations
        end subroutine
        
        !> Set activities of dissolved solid species with validation
        subroutine set_act_diss_solids(this,act_diss_solids)
        class(aqueous_chemistry_c) :: this                                 !< Aqueous chemistry object (polymorphic)
        real(kind=8), intent(in) :: act_diss_solids(:)                     !< Activities of dissolved solid species
        if (size(act_diss_solids)/=this%aq_phase%num_species-this%aq_phase%wat_flag) then
            error stop "Dimension error in activity of dissolved solids"   !< Dimension validation error
        end if
        this%activities(this%ind_var_act_species)=act_diss_solids          !< Assign activities (workaround implementation)
        end subroutine
        
        !> Set concentrations of secondary variable activity species across all phases
        subroutine set_conc_var_act_species(this,conc_var_act_species)
        !> variable activity species ordered in:
        !! aqueous primary species
        !! solid primary species (master surface site)
        !! aqueous secondary variable activity species
        !! surface complex secondaries (cation-exchanged)
        !! variable activity gases in equilibrium
        
        !! FALTA COMPROBAR CONCENTRACIONES NEGATIVAS                         !< TODO: Check for negative concentrations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_var_act_species(:)            !< Variable activity species concentrations
            if (size(conc_var_act_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species) then
                error stop "Dimension error in concentration of variable activity species" !< Dimension validation
            else                                                           !< If dimensions are correct
                !< Layout of conc_var_act_species (matches get_conc_var_act / get_conc_old_var_act_new_rz):
                !<   1..n_aq_prim                       : aqueous primary
                !<   n_aq_prim+1..n_prim                : solid primary (master surface site)
                !<   n_prim+1..n_prim+n_aq_sec          : aqueous secondary variable activity
                !<   n_prim+n_aq_sec+1..+num_exch_cats  : surface complex secondaries (cations)
                !<   <gas slots after that>
                associate(spc => this%solid_chemistry%reactive_zone%speciation_alg, &
                          rz  => this%solid_chemistry%reactive_zone, &
                          sc  => this%solid_chemistry)
                !< Aqueous primary
                this%concentrations(this%ind_var_act_species(1:spc%num_aq_prim_species))=&
                    conc_var_act_species(1:spc%num_aq_prim_species)
                !< Aqueous secondary variable activity
                if (spc%num_aq_sec_var_act_species>0) then
                    this%concentrations(this%ind_var_act_species(spc%num_aq_prim_species+1:spc%num_aq_var_act_species))=&
                        conc_var_act_species(spc%num_prim_species+1:spc%num_prim_species+spc%num_aq_sec_var_act_species)
                end if
                if (rz%cat_exch_zone%num_surf_compl>0) then
                    !< Solid primary (master/free surface site) at slot n_prim
                    sc%concentrations(sc%mineral_zone%num_minerals+1)=&
                        conc_var_act_species(spc%num_prim_species)
                    !< Surface complex secondaries (cation-exchanged)
                    sc%concentrations(sc%mineral_zone%num_minerals+2:sc%num_solids)=&
                        conc_var_act_species(spc%num_prim_species+spc%num_aq_sec_var_act_species+1:&
                            spc%num_prim_species+spc%num_aq_sec_var_act_species+rz%cat_exch_zone%num_exch_cats)
                end if
                if (associated(this%gas_chemistry)) then                   !< If gas chemistry is present
                    this%gas_chemistry%concentrations(rz%gas_phase%num_gases_eq_cst_act+1:&
                        rz%gas_phase%num_gases_eq)=&
                        conc_var_act_species(spc%num_aq_var_act_species+rz%cat_exch_zone%num_surf_compl+1:&
                            spc%num_var_act_species-rz%gas_phase%num_gases_kin_var_act)
                end if
                end associate
            end if
        end subroutine
        
        !> Set concentrations of secondary variable activity species
        subroutine set_conc_sec_var_act_species(this,c2v) !> sets concentration of secondary variable activity species
        !> c2v must be ordered in:
        !! 1. aqueous secondary variable activity species
        !! 2. surface complexes
        !! 3. variable activity gases in equilibrium
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c2v(:)                            !< Secondary variable activity concentrations
            if (size(c2v)/=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions) then
                error stop "Dimension error in concentration of secondary variable activity species"
            else if (size(c2v)==this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species) then
                this%concentrations(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species&
                    +1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))=c2v
            else
                this%concentrations(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species&
                    +1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))=&
                    c2v(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species)
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                    this%solid_chemistry%num_solids)=&
                c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
                if (associated(this%gas_chemistry)) then
                    this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=&
                        c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
                end if
            end if
        end subroutine
        
        !> \brief Sets concentration of secondary aqueous species
        !!
        !! \details Assigns the provided concentration array to the secondary
        !! aqueous species positions in the concentrations vector.
        !!
        !! @param c2aq concentration array of secondary aqueous species
        !!
        subroutine set_conc_sec_aq_species(this,c2aq) !> sets concentration of secondary aqueous species
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2aq(:) !> concentration of secondary aqueous species
            !if (present(c2aq)) then
                if (size(c2aq)/=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)then
                    error stop "Dimension error in concentration of secondary aqueous species"
                else
                    this%concentrations(this%ind_sec_species)=c2aq
                end if
        end subroutine
        
        
        !subroutine update_conc_sec_var_act_species(this,c2v) !> updates concentration of secondary variable activity species & 
        !    !! assuming there are aqueous and solid species
        
        !subroutine update_conc_sec_aq_var_act_species(this,c2v_aq)!> updates concentration of secondary aqueous variable activity &
        !    !! species
        
        !> \brief Sets concentrations of all secondary species
        !!
        !! \details Distributes the secondary species concentration array into
        !! aqueous, mineral, surface complexation, and gas sub-arrays.
        !!
        !! @param c2 concentration of secondary species ordered as: aqueous complexes,
        !!        var act minerals eq, cst act minerals eq, surface complexes,
        !!        var act gases eq, cst act gases eq
        !!
        subroutine set_conc_sec_species(this,c2)
            !> c2 must be ordered in:
            !! 1. aqueous complexes
            !! 2. var act minerals in equilibrium
            !! 3. cst act minerals in equilibrium
            !! 4. surface complexes
            !! 5. var act gases in equilibrium
            !! 6. cst act gases in equilibrium
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            real(kind=8), intent(in) :: c2(:) !> concentration of secondary species

            if (size(c2)/=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species) then
                error stop "Dimension error in 'set_conc_sec_species'"
            else if (this%solid_chemistry%reactive_zone%speciation_alg%flag_comp.eqv..false.) then
                this%concentrations(this%ind_sec_species)=c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals_kin+1:&
                    this%solid_chemistry%mineral_zone%num_minerals)=c2(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    this%solid_chemistry%reactive_zone%num_minerals)
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                        this%solid_chemistry%num_solids)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                        this%solid_chemistry%reactive_zone%num_minerals+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                        this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
                end if
                if (associated(this%gas_chemistry)) then
                    this%gas_chemistry%concentrations(1:this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species)
                    this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)
                end if
            end if
        end subroutine
       

        !> \brief Computes component concentrations from variable activity species
        !!
        !! \details Multiplies the component matrix by the variable activity species
        !! concentration vector to obtain the component concentrations.
        !!
        !! @param c_nc variable activity species concentrations
        !! @return conc_comp component concentrations
        !!
        function compute_conc_comp(this,c_nc) result(conc_comp)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_nc(:) !> variable activity species
            real(kind=8), allocatable :: conc_comp(:) !> component concentrations 
            conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_nc)
        end function
        
        
        !> \brief Computes component concentrations including constant activity species
        !!
        !! \details Multiplies the constant activity component matrix by the full
        !! species concentration vector to obtain component concentrations.
        !!
        !! @param conc concentration of all species
        !! @return conc_comp component concentrations
        !!
        function compute_conc_comp_cst_act(this,conc) result(conc_comp)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:) !> concentration of all species
            real(kind=8), allocatable :: conc_comp(:)
            conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act,conc)
        end function
        
      
        !> \brief Sets the ionic strength value
        !!
        !! @param ionic_strength ionic strength value [mol/kg]
        !!
        subroutine set_ionic_strength(this,ionic_strength)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: ionic_strength
            this%ionic_strength=ionic_strength
        end subroutine
        
        !> \brief Computes ionic strength of the aqueous solution
        !!
        !! \details Calculates the ionic strength (I) using the formula:
        !! \f[
        !! I = \frac{1}{2} \sum_{i=1}^{n_{species}} m_i z_i^2
        !! \f]
        !! where:
        !! - \f$ m_i \f$ is the molality of species i [mol/kg water]
        !! - \f$ z_i \f$ is the charge (valence) of species i [-]
        !! - The summation is over all aqueous species
        !!
        !! The ionic strength is a measure of the total concentration of ions in solution
        !! and is critical for calculating activity coefficients in non-ideal solutions.
        !! Higher ionic strength leads to greater deviations from ideal behavior.
        !!
        !! \note Concentrations must be in molalities (mol/kg water), not molarities (mol/L).
        !!       The ionic strength is stored in this%ionic_strength.
        !!
        !! \note this%concentrations follows THIS water's species order (the order of its water type
        !!       file), which can differ from the aqueous phase order: concentration index i belongs to
        !!       the aqueous species indices_aq_phase(i), whose valence must be used. Pairing
        !!       concentrations(i) with aq_species(i) would, e.g., weight the water concentration
        !!       (55.5 mol/kg) with the charge of an ion (I ~ 100 mol/kg, pH < 0).
        !!
        !! \warning This subroutine assumes all species in aq_phase%aq_species have valid
        !!          valence values. Species with zero valence (neutral species) contribute
        !!          zero to the ionic strength.
        !!
        !! \see compute_log_act_coeffs, set_ionic_strength
        !!
        subroutine compute_ionic_strength(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            integer(kind=4) :: i                                            !< Loop counter for aqueous species (concentration index)
            
            !> \subsection init_ionic_strength Initialize ionic strength to zero
            !! Start with zero and accumulate contributions from each species.
            this%ionic_strength=0d0                                         !< Initialize ionic strength to zero [mol/kg]
            
            !> \subsection loop_species Loop through all aqueous species
            !! Calculate contribution of each species to ionic strength: molality times squared valence
            do i=1,this%aq_phase%num_species                                !< Iterate over all aqueous species
                !> Add contribution of the species stored at concentration index i: molality times squared valence
                !! - this%concentrations(i): molality of the species at concentration index i [mol/kg water]
                !! - this%aq_phase%aq_species(this%indices_aq_phase(i))%valence: charge of THAT species [-]
                !! - Squared valence ensures positive contribution regardless of charge sign
                this%ionic_strength=this%ionic_strength+&
                    this%concentrations(i)*this%aq_phase%aq_species(this%indices_aq_phase(i))%valence**2
                                                                            !< Accumulate molality times squared valence term
            end do
            
            !> \subsection scale_factor Apply the 1/2 factor
            !! The ionic strength definition includes a factor of 0.5:
            !! I = (1/2) * sum(c_i * z_i^2)
            this%ionic_strength=0.5*this%ionic_strength                     !< Multiply by 1/2 to get final ionic strength [mol/kg]
        end subroutine
        
        !> \brief Computes this water's logarithmic equilibrium constants for all reactions
        !!
        !! \details This subroutine builds this water's per-reaction logÃ¢â€šÂÃ¢â€šâ‚¬(K(T)) at its own
        !! temperature for ALL reactions of the reactive zone (equilibrium block followed by the
        !! linear, redox and mineral kinetic blocks) and stores them in this%logK. For each
        !! reaction the temperature dependence follows the PHREEQC analytical expression:
        !! \f[
        !! \log_{10}(K(T)) = a_1 + a_2 T + \frac{a_3}{T} + a_4 \log_{10}(T) + \frac{a_5}{T^2} + a_6 T^2
        !! \f]
        !! where:
        !! - T is absolute temperature [K]
        !! - coeffs_logK_T = [a_1, a_2, a_3, a_4, a_5, a_6] are the 6 reaction coefficients
        !!
        !! When the 6 temperature coefficients are not available for a reaction, the value falls
        !! back to the reference-temperature constant logÃ¢â€šÂÃ¢â€šâ‚¬(eq_cst_ref) (0 when no equilibrium
        !! constant is defined).
        !!
        !! **Algorithm**:
        !! 1. Fill the equilibrium block of this%logK (reactions mapped via ind_eq_reacts)
        !! 2. Append the linear, redox and mineral kinetic blocks (chemical-system order)
        !!
        !! \note The temperature is taken from this%temp which must be set beforehand.
        !!       The computed logK values are per-water state; the shared reaction objects
        !!       keep only the reference-temperature eq_cst_ref.
        !!
        !! \see reaction_c::compute_logK_dep_T_PHREEQC, set_temp
        !!
        subroutine compute_logK_PHREEQC_aq_chem(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            integer(kind=4) :: i                                            !< Loop counter
            integer(kind=4) :: n_eq,n_lin,n_redox,n_min_kin,offset          !< logK block sizes and running offset
            
            !> Build this water's per-reaction log10(K_eq(T)) at its own temperature for ALL
            !! reactions, ordered like get_S_nc: the equilibrium reactions first (aqueous
            !! complexation + redox AND heterogeneous mineral/gas/exchange, mapped to the
            !! chemical system via ind_eq_reacts), followed by the kinetic reactions (linear,
            !! redox and mineral, in chemical-system order). The kinetic block is included so
            !! that saturation indices of kinetic minerals can use this per-water T-dependent
            !! constant; it is NOT used to build logK_tilde/logK_star (those only involve the
            !! equilibrium block). The shared reaction objects keep only the reference-T
            !! eq_cst_ref; the T-dependent value is per-water state.
            associate(rz => this%solid_chemistry%reactive_zone, &
                      cs => this%solid_chemistry%reactive_zone%chem_syst)
                n_eq      = rz%speciation_alg%num_eq_reactions
                n_lin     = cs%num_lin_kin_reacts
                n_redox   = cs%num_redox_kin_reacts
                n_min_kin = cs%num_minerals_kin
                if (allocated(this%logK)) deallocate(this%logK)
                allocate(this%logK(n_eq+n_lin+n_redox+n_min_kin))
                !> Equilibrium reactions
                do i=1,n_eq
                    this%logK(i)=logK_react_T(cs%eq_reacts(rz%ind_eq_reacts(i)),this%temp)
                end do
                offset=n_eq
                !> Linear kinetic reactions
                do i=1,n_lin
                    this%logK(offset+i)=logK_react_T(cs%lin_kin_reacts(i),this%temp)
                end do
                offset=offset+n_lin
                !> Redox kinetic reactions
                do i=1,n_redox
                    this%logK(offset+i)=logK_react_T(cs%redox_kin_reacts(i),this%temp)
                end do
                offset=offset+n_redox
                !> Mineral kinetic reactions (chemical-system order)
                do i=1,n_min_kin
                    this%logK(offset+i)=logK_react_T(cs%min_kin_reacts(i),this%temp)
                end do
            end associate
        end subroutine

            !> log10(K(T)) for one reaction at a water's temperature temp [K].
            !! At (approximately) the reference temperature there is NO temperature dependence to
            !! resolve, so it returns log10(eq_cst_ref) directly. Away from the reference temperature
            !! it uses, in order of preference: the PHREEQC analytic expression (6 coefficients
            !! available), the integrated Van't Hoff equation (delta_h_ref set), and otherwise
            !! log10(eq_cst_ref) (0 when no equilibrium constant is defined) with a one-time warning,
            !! since that silently leaves the reaction at 25 C. Module-level so that both the
            !! per-water logK blocks (compute_logK_PHREEQC_aq_chem) and the saturation indices of
            !! minerals outside the water's reactive zone (compute_saturation_eq_min) use it.
            function logK_react_T(rct,temp) result(lk)
                class(reaction_c), intent(in) :: rct !< reaction (any subclass)
                real(kind=8), intent(in) :: temp !< temperature of the water [K]
                real(kind=8) :: lk !< log10 of the equilibrium constant
                real(kind=8), parameter :: T_ref=298.15d0 !< reference temperature [K] (25 C)
                real(kind=8), parameter :: T_tol=1d-1 !< tolerance [K] to treat the water as being at the reference temperature
                if (abs(temp-T_ref)>=T_tol) then
                    if (allocated(rct%coeffs_logK_T)) then
                        if (size(rct%coeffs_logK_T)==6 .and. any(rct%coeffs_logK_T/=0d0)) then
                            lk=rct%compute_logK_dep_T_PHREEQC(temp)
                            return
                        end if
                    end if
                    if (rct%delta_h_ref/=0d0 .and. rct%eq_cst_ref>0d0) then
                        lk=rct%compute_logK_dep_T_Vant_Hoff(temp)
                        return
                    end if
                    if (rct%eq_cst_ref>0d0) call warn_no_T_data(rct,temp)
                end if
                if (rct%eq_cst_ref>0d0) then
                    lk=log10(rct%eq_cst_ref)
                else
                    lk=0d0
                end if
            end function

            !> Notify (once per reaction name) that an off-reference-temperature water uses the
            !! database log K unchanged because the database gives neither analytic coefficients
            !! nor delta_h (intended for entries already tabulated at the water's temperature).
            subroutine warn_no_T_data(rct,temp)
                class(reaction_c), intent(in) :: rct
                real(kind=8), intent(in) :: temp !< temperature of the water [K]
                integer(kind=4), parameter :: max_warned=64
                character(len=256), save :: warned(max_warned)
                integer(kind=4), save :: num_warned=0
                integer(kind=4) :: w
                do w=1,num_warned
                    if (warned(w)==rct%name) return
                end do
                if (num_warned<max_warned) then
                    num_warned=num_warned+1
                    warned(num_warned)=rct%name
                end if
                write(0,'(A,A,A,F7.2,A)') ' NOTE: reaction ''', trim(rct%name), &
                    ''' has no temperature data (analytic coefficients or delta_h): its database log K is used unchanged for a water at ', &
                    temp-273.15d0, ' C'
            end subroutine
       
        !> \brief Associates a solid chemistry object to this aqueous chemistry
        !!
        !! \details This subroutine creates a pointer association between the aqueous
        !! chemistry object and a solid chemistry object (containing minerals and their
        !! reactive zones). This coupling is necessary for:
        !! - Mineral dissolution/precipitation reactions
        !! - Solid-aqueous equilibria
        !! - Reactive transport with heterogeneous reactions
        !! - Updating equilibrium constants for mineral reactions
        !!
        !! **Algorithm**:
        !! 1. Verify that solid_chemistry%reactive_zone is associated
        !! 2. If valid, create pointer association this%solid_chemistry => solid_chemistry
        !! 3. If invalid, terminate with error message
        !!
        !! \param[in] solid_chemistry The solid chemistry object to associate (target attribute required)
        !!
        !! \note The solid_chemistry argument must have the target attribute so it can be
        !!       pointed to. The reactive_zone component must already be associated before
        !!       calling this subroutine.
        !!
        !! \warning This subroutine will terminate the program if the solid chemistry
        !!          object is not properly initialized with a reactive zone.
        !!
        !! \see set_gas_chemistry, solid_chemistry_c, reactive_zone_c
        !!
        subroutine set_solid_chemistry(this,solid_chemistry)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object to modify
            class(solid_chemistry_c), intent(in), target :: solid_chemistry !< [in] Solid chemistry object to associate (must be target)
            logical :: new_target                                           !< .true. if the water is pointed to a DIFFERENT solid chemistry object
            
            !> \subsection validate_reactive_zone Validate solid chemistry initialization
            !! Check that the solid chemistry object has been properly initialized with
            !! a reactive zone before creating the association. The reactive zone contains
            !! information about which minerals are present and their reaction kinetics.
            if (associated(solid_chemistry%reactive_zone)) then            !< Check if reactive zone is associated
                new_target=.not. associated(this%solid_chemistry,solid_chemistry)
                !> Create pointer association to solid chemistry
                !! This allows the aqueous chemistry object to access mineral data,
                !! equilibrium reactions involving minerals, and kinetic rate laws
                this%solid_chemistry=>solid_chemistry                      !< Point to the solid chemistry object
                !> The kinetic reaction set of a water is defined by its solid chemistry (the kinetic
                !> minerals of its mineral zone), so pointing the water to a different solid chemistry
                !> object changes it: rebuild indices_rk, and indices_rk_old with it (both sets carry the
                !> kinetic minerals of the CURRENT solid chemistry, see get_rk_old). Only for waters whose
                !> kinetic set was already built once (set_indices_rk at initialisation allocates
                !> indices_rk%cols, possibly with zero columns): while a water (type) is being read its
                !> solid chemistry is not set up yet. Nested ifs: Fortran does not short-circuit, and
                !> mineral_zone has no default initialisation on a fresh solid chemistry.
                if (new_target) then
                    if (allocated(this%indices_rk%cols) .and. associated(this%aq_phase)) then
                        if (associated(solid_chemistry%mineral_zone) .and. associated(solid_chemistry%reactive_zone%chem_syst)) then
                            call this%set_indices_rk()
                            this%indices_rk_old=this%indices_rk
                        end if
                    end if
                end if
            else
                !> Terminate if reactive zone is not associated
                !! This is a critical error because mineral reactions cannot be computed
                !! without reactive zone information
                error stop "Solid chemistry object is not associated to a reactive zone"
                                                                            !< Fatal error: invalid solid chemistry
            end if
        end subroutine

        !> \brief Associates a solid chemistry object representing the previous time step
        !!
        !! \details Creates a pointer association to a solid chemistry object
        !! representing the state at the previous time step. If not provided,
        !! defaults to the current solid chemistry, and the previous-step index arrays
        !! (ind_prim_species_old, ind_sec_species_old: the OLD zone's layout, read e.g. by
        !! get_conc_old_old_prim_new_rz) are seeded from the current ones, so every water
        !! (the external ones included, which the per-step update_attributes never visits) has a
        !! consistent previous-step state from the start.
        !!
        !! @param solid_chemistry_old optional solid chemistry object from previous time step;
        !!        defaults to current solid chemistry. When given, the caller is responsible for
        !!        the matching old index arrays (see copy_aq_chem_Lagr).
        !!
        subroutine set_solid_chemistry_old(this,solid_chemistry_old)
            implicit none
            class(aqueous_chemistry_c) :: this
            class(solid_chemistry_c), intent(in), target, optional :: solid_chemistry_old

            if (present(solid_chemistry_old)) then
                if (associated(solid_chemistry_old%reactive_zone)) then
                    this%solid_chemistry_old=>solid_chemistry_old
                else
                    error stop "solid chemistry object is not associated to a reactive zone"
                end if
            else
                this%solid_chemistry_old=>this%solid_chemistry
                if (allocated(this%ind_prim_species)) this%ind_prim_species_old=this%ind_prim_species
                if (allocated(this%ind_sec_species)) this%ind_sec_species_old=this%ind_sec_species
                this%indices_rk_old=this%indices_rk
            end if
        end subroutine
        
        !> \brief Associates a gas chemistry object to this aqueous chemistry
        !!
        !! \details This subroutine creates a pointer association between the aqueous
        !! chemistry object and a gas chemistry object (containing gases and their
        !! reactive zones). This coupling is necessary for:
        !! - Gas dissolution/exsolution reactions (e.g., COÃ¢â€šâ€š, OÃ¢â€šâ€š, CHÃ¢â€šâ€ž)
        !! - Gas-aqueous equilibria (Henry's law)
        !! - Reactive transport with gas phase
        !! - Computing partial pressures and fugacities
        !!
        !! **Algorithm**:
        !! 1. Verify that gas_chemistry%reactive_zone is associated
        !! 2. If valid, create pointer association this%gas_chemistry => gas_chemistry
        !! 3. If invalid, terminate with error message
        !!
        !! \param[in] gas_chemistry The gas chemistry object to associate (target attribute required)
        !!
        !! \note The gas_chemistry argument must have the target attribute so it can be
        !!       pointed to. The reactive_zone component must already be associated before
        !!       calling this subroutine.
        !!
        !! \warning This subroutine will terminate the program if the gas chemistry
        !!          object is not properly initialized with a reactive zone.
        !!
        !! \see set_solid_chemistry, gas_chemistry_c, reactive_zone_c
        !!
        subroutine set_gas_chemistry(this,gas_chemistry)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object to modify
            class(gas_chemistry_c), intent(in), target :: gas_chemistry    !< [in] Gas chemistry object to associate (must be target)
            
            !> \subsection validate_gas_reactive_zone Validate gas chemistry initialization
            !! Check that the gas chemistry object has been properly initialized with
            !! a reactive zone before creating the association. The reactive zone contains
            !! information about which gases are present and their Henry's law constants.
            if (associated(gas_chemistry%reactive_zone)) then              !< Check if reactive zone is associated
                !> Create pointer association to gas chemistry
                !! This allows the aqueous chemistry object to access gas phase data,
                !! equilibrium reactions involving gases (Henry's law), and gas partial pressures
                this%gas_chemistry=>gas_chemistry                          !< Point to the gas chemistry object
            else
                !> Terminate if reactive zone is not associated
                !! This is a critical error because gas-aqueous equilibria cannot be computed
                !! without reactive zone information
                error stop "gas chemistry object is not associated to a reactive zone"
                                                                            !< Fatal error: invalid gas chemistry
            end if
        end subroutine
        
     
        !> \brief Sets the pH value
        !!
        !! \details Assigns the pH to the aqueous chemistry object.
        !! Defaults to 7.0 if not provided.
        !!
        !! @param pH optional pH value; defaults to 7.0
        !!
        subroutine set_pH(this,pH)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: pH
            if (present(pH)) then
                this%pH=pH
            else
                this%pH=7d0 !> default
            end if
        end subroutine
        
        !> \brief Computes pH from the activity of the proton species
        !!
        !! \details Calculates pH = -log10(a(H+)) using the activity of the
        !! proton species in the aqueous phase.
        !!
        subroutine compute_pH(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (this%aq_phase%ind_prot/=0) then !> chapuza
                this%pH=-log10(this%activities(this%indices_aq_species(this%aq_phase%ind_prot)))
            end if
        end subroutine
        
        !> \brief Sets the electron activity (pe) value
        !!
        !! \details Assigns the pe value to the aqueous chemistry object.
        !! Defaults to 4.0 if not provided.
        !!
        !! @param pe optional pe value; defaults to 4.0
        !!
        subroutine set_pe(this,pe)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: pe
            if (present(pe)) then
                this%pe=pe
            else
                this%pe=4d0 !> default
            end if
        end subroutine
        
        !> \brief Sets the concentration of a single aqueous species by index
        !!
        !! \details Assigns a concentration value to a specific aqueous species
        !! identified by its index in the aqueous phase.
        !!
        !! @param conc_sp concentration value of the species
        !! @param sp_ind index of the species in the aqueous phase
        !!
        subroutine set_conc_single_species(this,conc_sp,sp_ind)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_sp !> concentration of species
            integer(kind=4), intent(in) :: sp_ind !> index of species
            if (sp_ind<1 .or. sp_ind>this%aq_phase%num_species) then
                error stop "Error in aqueous species index"
            else if (.not. allocated(this%concentrations)) then
                error stop "Aqueous species concentrations must be allocated"
            else if (conc_sp<0d0) then
                error stop "Concentrations cannot be negative"
            end if
            this%concentrations(sp_ind)=conc_sp
        end subroutine
        
      
        !> \brief Allocates arrays for kinetic and equilibrium reaction rates
        !!
        !! \details Allocates Rk, rk_old, rk_new, rk_mean, cum_rxn_amnt_kin,
        !! re_mean, and Re arrays for aqueous, solid, and gas chemistry objects.
        !! All arrays are initialised to zero.
        !!
        subroutine allocate_rxn_rates_and_amounts(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (.not. associated(this%solid_chemistry)) then
                error stop "solid_chemistry pointer not associated in allocate_rxn_rates_and_amounts"
            end if
            if (.not. associated(this%solid_chemistry%reactive_zone)) then
                error stop "reactive_zone pointer not associated in allocate_rxn_rates_and_amounts"
            end if
            if (associated(this%solid_chemistry%reactive_zone%chem_syst)) then
                associate(n_kin => this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                                   this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts)
                allocate(this%rxn_amnt_kin(n_kin))
                allocate(this%rxn_rate_kin_old(n_kin))
                allocate(this%rk_old_old(n_kin))
                allocate(this%rk_old_old_old(n_kin))
                allocate(this%rxn_rate_kin_avg(n_kin))
                allocate(this%rxn_rate_kin(n_kin))
                allocate(this%cum_rxn_amnt_kin(n_kin))
                !allocate(this%contrib_rxn_amnt_kin(n_kin))
                end associate
                this%rxn_amnt_kin=0d0 !> by default
                this%rxn_rate_kin_old=0d0 !> by default
                this%rxn_rate_kin=0d0 !> by default
                this%rk_old_old=0d0 !> by default
                this%rk_old_old_old=0d0 !> by default
                this%rxn_rate_kin_avg=0d0 !> by default
                this%cum_rxn_amnt_kin=0d0 !> by default
                !this%contrib_rxn_amnt_kin=0d0 !> by default
            else
                error stop "Chemical system not associated to reactive zone"
            end if
            allocate(this%rxn_rate_eq(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%rxn_rate_eq_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%rxn_rate_eq_old_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%rxn_rate_eq_old_old_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%rxn_amnt_eq(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            !allocate(this%contrib_rxn_amnt_eq(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%cum_rxn_amnt_eq(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            this%rxn_rate_eq=0d0 !> by default
            this%rxn_rate_eq_old=0d0 !> by default
            this%rxn_rate_eq_old_old=0d0 !> by default
            this%rxn_rate_eq_old_old_old=0d0 !> by default
            this%rxn_amnt_eq=0d0 !> by default
            !this%contrib_rxn_amnt_eq=0d0 !> by default
            this%cum_rxn_amnt_eq=0d0 !> by default
            if (associated(this%solid_chemistry)) then
                if (.not. allocated(this%solid_chemistry%rxn_rate_eq)) then
                    allocate(this%solid_chemistry%rxn_rate_eq(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%rxn_rate_eq_old(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%rxn_rate_eq_old_old(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%rxn_rate_eq_old_old_old(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%rxn_amnt_eq(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    ! allocate(this%solid_chemistry%contrib_rxn_amnt_eq(this%solid_chemistry%reactive_zone%num_minerals+&
                    !     this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%cum_rxn_amnt_eq(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    this%solid_chemistry%rxn_rate_eq=0d0 !> by default
                    this%solid_chemistry%rxn_rate_eq_old=0d0 !> by default
                    this%solid_chemistry%rxn_rate_eq_old_old=0d0 !> by default
                    this%solid_chemistry%rxn_rate_eq_old_old_old=0d0 !> by default
                    this%solid_chemistry%rxn_amnt_eq=0d0 !> by default
                    !this%solid_chemistry%contrib_rxn_amnt_eq=0d0 !> by default
                    this%solid_chemistry%cum_rxn_amnt_eq=0d0 !> by default
                end if
                if (.not. allocated(this%solid_chemistry%rxn_amnt_kin)) then
                    associate(n_min_kin => this%solid_chemistry%mineral_zone%num_minerals_kin)
                    allocate(this%solid_chemistry%rxn_amnt_kin(n_min_kin))
                    allocate(this%solid_chemistry%rxn_rate_kin_old(n_min_kin))
                    allocate(this%solid_chemistry%rk_old_old(n_min_kin))
                    allocate(this%solid_chemistry%rk_old_old_old(n_min_kin))
                    allocate(this%solid_chemistry%rxn_rate_kin_avg(n_min_kin))
                    allocate(this%solid_chemistry%rxn_rate_kin(n_min_kin))
                    allocate(this%solid_chemistry%cum_rxn_amnt_kin(n_min_kin))
                    !allocate(this%solid_chemistry%contrib_rxn_amnt_kin(n_min_kin))
                    end associate
                    this%solid_chemistry%rxn_rate_kin=0d0 !> by default
                    this%solid_chemistry%rxn_rate_kin_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old_old=0d0 !> by default
                    this%solid_chemistry%rxn_rate_kin_avg=0d0 !> by default
                    this%solid_chemistry%rxn_amnt_kin=0d0 !> by default
                    this%solid_chemistry%cum_rxn_amnt_kin=0d0 !> by default
                    !this%solid_chemistry%contrib_rxn_amnt_kin=0d0 !> by default
                end if
            else
                error stop "Solid chemistry not associated to aqueous chemistry"
            end if
            if (associated(this%gas_chemistry)) then
                if (.not. allocated(this%gas_chemistry%rxn_rate_eq)) then
                    allocate(this%gas_chemistry%rxn_rate_eq(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%rxn_rate_eq_old(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%rxn_rate_eq_old_old(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%rxn_rate_eq_old_old_old(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%rxn_amnt_eq(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    !allocate(this%gas_chemistry%contrib_rxn_amnt_eq(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    this%gas_chemistry%rxn_rate_eq=0d0 !> by default
                    this%gas_chemistry%rxn_rate_eq_old=0d0 !> by default
                    this%gas_chemistry%rxn_rate_eq_old_old=0d0 !> by default
                    this%gas_chemistry%rxn_rate_eq_old_old_old=0d0 !> by default
                    this%gas_chemistry%rxn_amnt_eq=0d0 !> by default
                    !this%gas_chemistry%contrib_rxn_amnt_eq=0d0 !> by default
                end if
            end if
        end subroutine
        
        !> \brief Associates an aqueous phase object to this aqueous chemistry
        !!
        !! \details This subroutine creates a pointer association between the aqueous
        !! chemistry object and an aqueous phase object. The aqueous phase contains
        !! information about:
        !! - Water properties (density, viscosity, dielectric constant)
        !! - Aqueous species present in the phase
        !! - Activity coefficient models
        !! - Temperature and pressure conditions
        !!
        !! **Algorithm**:
        !! 1. Create pointer association this%aq_phase => aq_phase
        !!
        !! \param[in] aq_phase The aqueous phase object to associate (target attribute required)
        !!
        !! \note The aq_phase argument must have the target attribute so it can be
        !!       pointed to. Unlike set_solid_chemistry and set_gas_chemistry, this
        !!       subroutine does not validate association status because aqueous phase
        !!       is always required for aqueous chemistry.
        !!
        !! \see set_solid_chemistry, set_gas_chemistry, aq_phase_c
        !!
        subroutine set_aq_phase(this,aq_phase)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object to modify
            class(aq_phase_c), intent(in), target :: aq_phase              !< [in] Aqueous phase object to associate (must be target)
            
            !> Create pointer association to aqueous phase
            !! This direct assignment establishes the connection without validation
            !! because aqueous chemistry always requires an aqueous phase
            this%aq_phase=>aq_phase                                        !< Point to the aqueous phase object
        end subroutine
        
      
        !> \brief Computes activity of an aqueous species from the mass action law
        !!
        !! \details Determines the activity of a specified aqueous species by
        !! solving the mass action law of the given equilibrium reaction.
        !!
        !! @param ind_act_aq_phase index of the species with unknown activity in aqueous phase
        !! @param ind_eq_rct index of the equilibrium reaction in the chemical system
        !! @param act_non_aq_sp optional activity of non-aqueous species (partial pressure if gas)
        !!
        subroutine compute_act_from_MAL(this,ind_act_aq_phase,ind_eq_rct,act_non_aq_sp) !> computes activity aqueous species from mass action law
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: ind_act_aq_phase !> index species unknown activity in aqueous phase
            !class(eq_reaction_c), intent(in) :: eq_react !> equilibrium reaction
            integer(kind=4), intent(in) :: ind_eq_rct !> index of equilibrium reaction in chemical system
            real(kind=8), intent(in), optional :: act_non_aq_sp !> activity of non-aqueous species (partial pressure if gas)
            
            real(kind=8) :: log_a,dot_prod
            integer(kind=4) :: ind_act_react,i
            integer(kind=8), allocatable :: react_indices(:)

            react_indices=this%aq_phase%get_indices_reaction(&
                this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct))

            dot_prod=0d0
            do i=1,size(react_indices)
                if (ind_act_aq_phase==react_indices(i)) then
                    ind_act_react=i
                else
                    dot_prod=dot_prod+this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%stoichiometry(i)*log10(this%activities(react_indices(i)))
                end if
            end do
            if (present(act_non_aq_sp)) then
                if (size(react_indices)==this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%num_species-1) then
                    dot_prod=dot_prod+log10(act_non_aq_sp)
                else
                    error stop "Dimension error in compute_act_from_MAL"
                end if
            end if
            log_a=(log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%eq_cst_ref)-dot_prod)/this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%stoichiometry(ind_act_react)
            this%activities(ind_act_aq_phase)=10**log_a
        end subroutine
        
        !> \brief Converts concentrations from molalities to molarities
        !!
        !! \details Applies density and salinity corrections to convert
        !! dissolved solids from molality to molarity units.
        !!
        subroutine compute_molarities(this) !> from molalities
            implicit none
            class(aqueous_chemistry_c) :: this
            
            ! integer(kind=4) :: i
            
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))*this%density*(1d0-this%salinity)
            
            !> Apply density correction to water concentration if water flag is set
            if (this%aq_phase%wat_flag==1) then                            !< Check if water is explicitly tracked as a species
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=this%concentrations(&
                    this%indices_aq_species(this%aq_phase%ind_wat))*this%density  !< Multiply water concentration by fluid density
            end if

        end subroutine
        
        !> Compute salinity from dissolved solids concentrations
        !! Mathematical Context: Salinity = TDS/(1+TDS) where TDS = ÃŽÂ£(ci Ãƒâ€” MWi) [kg_solute/kg_water]
        !! This converts molality-based concentrations to mass fraction salinity
        !> \brief Computes salinity (total dissolved solids mass fraction)
        !!
        !! \details This subroutine calculates the salinity of the aqueous solution as
        !! the mass fraction of dissolved solids. Salinity is defined as:
        !! \f[
        !! S = \frac{m_{solutes}}{m_{solutes} + m_{water}} = \frac{TDS}{1 + TDS}
        !! \f]
        !! where TDS (Total Dissolved Solids) is computed from:
        !! \f[
        !! TDS = \sum_{i=1}^{n_{diss}} m_i \cdot MW_i
        !! \f]
        !! where:
        !! - S = salinity [-] (mass fraction, dimensionless)
        !! - TDS = total dissolved solids [kg_solute/kg_water]
        !! - \f$ m_i \f$ = molality of dissolved species i [mol/kg_water]
        !! - \f$ MW_i \f$ = molecular weight of species i [kg/mol]
        !! - \f$ n_{diss} \f$ = number of dissolved solid species (excluding water)
        !!
        !! **Algorithm**:
        !! 1. Initialize TDS = 0
        !! 2. Loop over all dissolved species and accumulate mass: TDS += m_i * MW_i
        !! 3. Compute total fluid mass = 1 kg water + TDS kg solutes
        !! 4. Compute salinity = TDS / total_fluid_mass
        !!
        !! \note This subroutine assumes concentrations are in molalities [mol/kg_water].
        !!       Water itself is excluded from the sum (water is the solvent, not solute).
        !!
        !! \warning The calculation assumes 1 kg of water as the reference basis for molalities.
        !!
        !! \see compute_ionic_strength, compute_TDS
        !!
        subroutine compute_salinity(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            real(kind=8) :: TDS                                            !< Total Dissolved Solids [kg_solute/kg_water]
            real(kind=8) :: fluid_mass                                     !< Total fluid mass [kg]
            integer(kind=4) :: i                                           !< Loop index for dissolved solid species
            
            !> \subsection init_tds Initialize TDS accumulator
            !! Start with zero mass of dissolved solids
            TDS=0d0                                                         !< Initialize TDS to zero [kg_solute/kg_water]
            
            !> \subsection accumulate_mass Sum mass contributions from all dissolved species
            !! Loop over all dissolved solid species (ions, complexes) and add their
            !! mass contributions: mass_i = molality_i * molecular_weight_i
            do i=1,this%aq_phase%num_species-this%aq_phase%wat_flag         !< Loop over dissolved solids (water excluded)
                !> Add mass contribution for species i
                !! Extract concentration from full array using index mapping,
                !! then multiply by molecular weight to convert mol/kg -> kg/kg
                TDS=TDS+this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids(i))) * &
                    this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i))%molecular_weight
                                                                            !< TDS += m_i * MW_i [kg_solute/kg_water]
            end do
            
            !> \subsection compute_total_mass Calculate total fluid mass
            !! Total mass = mass of water (1 kg basis) + mass of all dissolved solutes
            fluid_mass=1d0+TDS                                             !< Total fluid mass = 1 kg water + TDS [kg]
            
            !> \subsection compute_salinity_value Compute salinity as mass fraction
            !! Salinity is the dimensionless ratio of solute mass to total fluid mass
            this%salinity=TDS/fluid_mass                                   !< S = TDS/(1+TDS) [-]
        end subroutine
        
        !> Compute alkalinity as charge balance contribution from all aqueous species
        !! Mathematical Context: Alkalinity = ÃŽÂ£(ci Ãƒâ€” alk_contrib_i) where alk_contrib is species charge contribution
        !! Alkalinity represents the acid-neutralizing capacity of the solution
        subroutine compute_alkalinity(this)                                 !< Subroutine assumes concentrations are in molalities
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(inout) :: this              !< Aqueous chemistry object (polymorphic, modified)
            real(kind=8) :: alkalinity                                     !< Local alkalinity accumulator variable
            integer(kind=4) :: i                                           !< Loop index for aqueous species
            alkalinity = 0d0                                               !< Initialize alkalinity to zero
            !> Sum alkalinity contributions from all aqueous species
            do i=1, this%aq_phase%num_species                              !< Loop over all aqueous species
                alkalinity = alkalinity + this%concentrations(this%indices_aq_species(i)) * this%aq_phase%aq_species(i)%alk_contrib
                                                                            !< Add contribution: concentration Ãƒâ€” alkalinity_factor
            end do
            this%alkalinity = alkalinity                                   !< Store computed alkalinity in object
        end subroutine
        
        !> Convert molarities to molalities using density and salinity corrections
        !! Mathematical Context: molality = molarity / (density * (1 - salinity))
        !! This accounts for volume changes due to dissolved solids and temperature effects
        subroutine compute_molalities(this)                                 !< Convert from molarities to molalities
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            
            ! real(kind=8) :: log_a,dot_prod                                !< Commented legacy variables
            ! integer(kind=4) :: ind_act_react,i                           !< Commented legacy variables
            
            !> Convert dissolved solids concentrations from molarity to molality
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))/(this%density*(1d0-this%salinity))
                                                                            !< molality = molarity / (density Ãƒâ€” (1-salinity))
            !> Convert water concentration using density correction only
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))/this%density
                                                                            !< water_molality = water_molarity / density
        end subroutine

        !> Reactive surfaces of the minerals of this water's solid chemistry on a per-mass-of-water
        !! (molality) basis, returned in react_surfaces; the solid chemistry itself is NOT modified
        !! (it is shared by all the waters pointing to it and keeps the input basis).
        !! Mathematical Context: the stored reactive surfaces are per litre of water, [m2 mineral/L water];
        !! a litre of solution of density rho [kg/L] and salinity S [kg solute/kg solution] holds
        !! rho*(1-S) kg of water, so, exactly as compute_molalities does for the solutes,
        !!   react_surf [m2/kg H2O] = react_surf [m2/L water] / (density * (1 - salinity)),
        !! with THIS water's density and salinity. The kinetic rate laws fed with it then yield rates
        !! per kg of water, consistent with molalities. Zero-sized result if the solid chemistry has no
        !! reactive surfaces.
        subroutine compute_react_surfaces_per_kg_water(this,react_surfaces) !< m2/L water -> m2/kg H2O
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (polymorphic)
            real(kind=8), allocatable, intent(out) :: react_surfaces(:)    !< Reactive surfaces per kg of water [m2 mineral/kg H2O], one per mineral of the mineral zone (solid_chemistry%react_surfaces order)
            if (associated(this%solid_chemistry)) then
                if (allocated(this%solid_chemistry%react_surfaces)) then
                    react_surfaces=this%solid_chemistry%react_surfaces/(this%density*(1d0-this%salinity))
                                                                            !< m2/kg H2O = (m2/L water) / (density * (1-salinity))
                    return
                end if
            end if
            allocate(react_surfaces(0))
        end subroutine

        !> Inverse of compute_react_surfaces_per_kg_water: reactive surfaces given per kg of water
        !! [m2 mineral/kg H2O] converted to per litre of water [m2 mineral/L water],
        !! react_surf * density * (1 - salinity) with THIS water's density and salinity (as
        !! compute_molarities does for the solutes), and STORED in the solid chemistry's react_surfaces
        !! (the input basis of that array).
        subroutine compute_react_surfaces_per_L_water(this,react_surfaces_per_kg) !< m2/kg H2O -> m2/L water
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: react_surfaces_per_kg(:)           !< Reactive surfaces per kg of water [m2 mineral/kg H2O]
            this%solid_chemistry%react_surfaces=react_surfaces_per_kg*this%density*(1d0-this%salinity)
                                                                            !< m2/L water = (m2/kg H2O) * density * (1-salinity)
        end subroutine

        !> Density of pure water at 1 atm as a function of temperature, stored in this%density [kg/L].
        !! Mathematical Context: Kell (1975) correlation, J. Chem. Eng. Data 20, 97-105, with t the
        !! temperature in Celsius (valid 0-150 C, max deviation ~1e-5 kg/L):
        !!   rho [kg/m3] = (999.83952 + 16.945176 t - 7.9870401e-3 t^2 - 46.170461e-6 t^3
        !!                  + 105.56302e-9 t^4 - 280.54253e-12 t^5) / (1 + 16.879850e-3 t)
        !! (rho(25 C) = 997.05 kg/m3 = 0.99705 kg/L, the value set_density uses as default). This is the
        !! density of the solvent; the molality/molarity conversions (compute_molalities,
        !! compute_molarities, compute_react_surfaces_per_*) combine it with the salinity to account
        !! for the dissolved solids.
        !! \param temp [optional] Temperature [K]; this%temp is used when absent (this%temp is NOT
        !!        modified when given: set it with set_temp). Outside the validity range the correlation
        !! is evaluated at the nearest bound (0 or 150 C) and a note is written once.
        subroutine compute_density_water(this,temp)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in), optional :: temp                     !< Temperature [K] (default: this%temp)
            real(kind=8), parameter :: T0=273.15d0                          !< 0 C in Kelvin [K]
            real(kind=8), parameter :: t_min=0d0, t_max=150d0               !< Validity range of the correlation [C]
            real(kind=8), parameter :: a0=999.83952d0, a1=16.945176d0, a2=-7.9870401d-3, a3=-46.170461d-6, &
                a4=105.56302d-9, a5=-280.54253d-12, b1=16.879850d-3        !< Kell (1975) coefficients (t in C, rho in kg/m3)
            logical, save :: range_warned=.false.                           !< one-time note for temperatures outside the validity range
            real(kind=8) :: t                                               !< Temperature [C]
            if (present(temp)) then
                t=temp-T0
            else
                t=this%temp-T0
            end if
            if (t<t_min .or. t>t_max) then
                if (.not. range_warned) then
                    write(0,'(A,F8.2,A)') ' NOTE: compute_density_water: temperature ', t, &
                        ' C outside the 0-150 C validity range of the Kell (1975) water density correlation; evaluated at the nearest bound'
                    range_warned=.true.
                end if
                t=min(max(t,t_min),t_max)
            end if
            this%density=(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*a5)))))/(1d0+b1*t)/1d3 !< kg/m3 -> kg/L
        end subroutine

        !> This water's Debye-Hückel parameters A and B at its temperature (this%temp) and pure-water
        !! density (this%density, set by compute_density_water), stored in this%params_aq_sol
        !! (params_aq_sol_s%compute_params_aq_sol: Bradley & Pitzer (1979) dielectric constant, PHREEQC
        !! convention). Like the temperature and the density, A and B are a property of the water's
        !! cell: call it after set_temp/compute_density_water, before the activity coefficients
        !! (compute_log_act_coeffs_aq_chem, the Picard speciations) are evaluated.
        subroutine compute_params_aq_sol_aq_chem(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            call this%params_aq_sol%compute_params_aq_sol(this%temp,this%density)
        end subroutine
        
        !> Get primary species activities for speciation calculations
        !! Mathematical Context: Primary activities are basis for computing secondary species via mass action law
        !! These are the independent variables in the chemical equilibrium system
        function get_a1(this) result(a1)                                   !< Function returns primary activities array
        implicit none                                                       !< No implicit variable declarations
        class(aqueous_chemistry_c), intent(in) :: this                     !< Aqueous chemistry object (input only)
        real(kind=8), allocatable :: a1(:)                                 !< Primary activities output array

        integer(kind=4) :: i                                               !< Loop index for primary species

        !> Allocate array for all primary species (aqueous + solid + gas)
        allocate(a1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
        !> Extract aqueous primary species activities
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            a1(i)=this%activities(this%ind_prim_species(i))                !< Copy activity from main activities array
        end do
        !> faltan los solidos y gases                                      !< TODO: Missing solid and gas primary species
        
        end function
        
        !> Get primary species concentrations for speciation calculations
        !! Mathematical Context: Primary concentrations are the basis for mass balance in speciation
        !! These must be conserved while secondary concentrations are computed from equilibrium
        function get_c1(this) result(c1)                                   !< Function returns primary concentrations array
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), allocatable :: c1(:)                             !< Primary concentrations output array
            
            integer(kind=4) :: i                                           !< Loop index for primary species
            
            !> Allocate array for all primary species (aqueous + solid + gas)
            allocate(c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            
            c1=1d0                                                          !< Initialize to 1.0 (temporary hack - "chapuza")
            !> Extract aqueous primary species concentrations
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                c1(i)=this%concentrations(this%ind_prim_species(i))     !< Copy concentration from main array
            end do
            !> Add surface complexation primary species if present
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                c1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+1)      !< Add surface site concentration
            end if
            !> Add gas phase primary species if present
            if  (associated(this%gas_chemistry)) then                      !< Check if gas chemistry is associated
                c1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=this%gas_chemistry%concentrations(&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
                    this%gas_chemistry%reactive_zone%gas_phase%num_species)/this%volume
                                                                            !< Convert gas concentrations to molar basis using volume
            end if
        end function
        
        !> Get primary species concentrations only (subset of get_c1)
        !! This function extracts only the aqueous portion of primary species for aqueous-only calculations
        function get_conc_old_old_prim_new_rz(this) result(c1_old_old)                                   !< Function returns primary concentrations array
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), allocatable :: c1_old_old(:)                             !< Primary concentrations output array
            
            integer(kind=4) :: i                                           !< Loop index for primary species
            
            !> Allocate array for all primary species (aqueous + solid + gas)
            allocate(c1_old_old(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            
            c1_old_old=1d0                                                          !< Initialize to 1.0 (temporary hack - "chapuza")
            !> Extract aqueous primary species concentrations
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                c1_old_old(i)=this%conc_old_old(this%ind_prim_species(i))     !< Copy concentration from main array
            end do
            !> Add surface complexation primary species if present
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                c1_old_old(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=&
                    this%solid_chemistry%conc_old_old(&
                    this%solid_chemistry%mineral_zone%num_minerals+1)      !< Add surface site concentration
            end if
            !> Add gas phase primary species if present
            if  (associated(this%gas_chemistry)) then                      !< Check if gas chemistry is associated
                c1_old_old(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                    this%gas_chemistry%conc_old_old(&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
                    this%gas_chemistry%reactive_zone%gas_phase%num_species)/this%volume
                                                                            !< Convert gas concentrations to molar basis using volume
            end if
        end function
        
        !> Get all secondary species concentrations in standardized order
        !! Mathematical Context: Secondary species concentrations computed from primary via mass action law
        !! Ordering: aq_sec_var_act, aq_sec_cst_act, minerals_eq_var_act, minerals_eq_cst_act, surf_complexes, gases_eq_var_act, gases_eq_cst_act
        function get_c2(this) result(c2)                                   !< Function returns secondary concentrations array
        !! ordered in:
        !! aqueous secondary var act species
        !! aqueous secondary cst act species
        !! var act minerals in equilibrium
        !! cst act minerals in equilibrium
        !! surface complexes
        !! var act gases in eq
        !! cst act gases in equilibrium
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), allocatable :: c2(:)                             !< Secondary concentrations output array
            integer(kind=4) :: i                                           !< Loop index for different species types
            
            !> Allocate array for all secondary species
            allocate(c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species))
        !> Aqueous secondary species (sorted in var act and cst act)
            c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)=this%concentrations(this%ind_sec_species)
                                                                            !< Copy aqueous secondary concentrations block
        !> Minerals in equilibrium (sorted in var act and cst act)
            do i=1,this%solid_chemistry%reactive_zone%num_minerals          !< Loop over equilibrium minerals
                c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)  !< Copy equilibrium mineral concentrations
            end do
        !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats  !< Loop over surface exchange categories
                c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    this%solid_chemistry%reactive_zone%num_minerals+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i)    !< Copy surface complex concentrations
            end do
        !> Gases in equilibrium
            if (associated(this%gas_chemistry)) then                       !< Check if gas chemistry is present
                !> Gases with var act
                do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act  !< Loop over variable activity gases
                    c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+i)=this%gas_chemistry%concentrations(&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)  !< Copy variable activity gas concentrations
                end do
                !> Gases with cst act
                do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act  !< Loop over constant activity gases
                    c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)=this%gas_chemistry%concentrations(i)
                                                                            !< Copy constant activity gas concentrations
                end do
            end if
        end function
        
        !> Get secondary variable activity species concentrations only
        !! Mathematical Context: These are secondary species whose activities vary with solution composition
        !! Used in speciation calculations where activity coefficients must be computed iteratively
        function get_c2v(this) result(c2v)                               !< Function returns secondary variable activity concentrations
        !! ordered in:
        !! sec var act aq species
        !! var act minerals in eq
        !! surf complexes
        !! var act gases in eq
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), allocatable :: c2v(:)                           !< Secondary variable activity concentrations output
            
            integer(kind=4) :: i
            
            allocate(c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                c2v(i)=this%concentrations(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                this%solid_chemistry%reactive_zone%num_minerals_var_act)=this%solid_chemistry%concentrations(&
                this%solid_chemistry%mineral_zone%num_minerals_kin+1:this%solid_chemistry%mineral_zone%num_minerals_kin+&
                this%solid_chemistry%reactive_zone%num_minerals_var_act)
            if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=&
                    this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                    this%solid_chemistry%num_solids)
            end if
            if (associated(this%gas_chemistry)) then
                !if (associated(this%solid_chemistry)) then
                    c2v(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%num_solids:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)=this%gas_chemistry%concentrations(&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
            end if
        end function
        
        !> \brief Gets concentrations of variable activity species per unit volume of water
        !!
        !! \details Assembles and returns variable activity species concentrations.
        !! Ordered in: aq primary, solid primary, gas primary, aq secondary var act,
        !! surface complexes, gases eq var act.
        !!
        !! @return conc_nc variable activity species concentrations
        !!
        function get_conc_var_act(this) result(conc_nc) !> gets conc var act species per unit volume of water
        !> we assume all minerals have cst act
        !> Ordered in:
        !! aqueous primary species
        !! solid primary species
        !! gas primary species
        !! aqueous secondary var act species
        !! surface complexes
        !! gases eq var act
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc_nc(:)
            
            integer(kind=4) :: i
            
            allocate(conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                conc_nc(i)=this%concentrations(this%ind_var_act_species(i))
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
                !! this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=this%concentrations(&
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                        this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1)
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=this%solid_chemistry%concentrations(&
                        this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%num_solids)
                end if
            if (associated(this%gas_chemistry)) then
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act)=this%gas_chemistry%concentrations(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
            end if
        end function

        !> \brief Previous-step concentrations of the variable activity species, laid out in the CURRENT
        !!        (new) reactive zone of the water.
        !!
        !! \details Returns the previous-step (k) concentrations c_v^k of the variable activity
        !! species in the column order of the reactive zone the water sits in NOW
        !! (solid_chemistry%reactive_zone), i.e. the previous-step state expressed in the basis the
        !! current-step speciation works in. Same ordering as get_conc_var_act:
        !!   1. aqueous primary species
        !!   2. solid primary species (exchange free site, if cation exchange)
        !!   3. gas primary species
        !!   4. aqueous secondary variable activity species
        !!   5. surface complexes (adsorbed cations)
        !!   6. equilibrium gases with variable activity
        !! Minerals are assumed to have constant activity and do not appear.
        !!
        !! The aqueous blocks are read from this%conc_old through the CURRENT index arrays
        !! (ind_prim_species, ind_sec_species): conc_old is stored in the water's own aqueous-species
        !! order, independent of the zone, so reading it through the current zone's indices re-expresses
        !! it in the new layout (the aqueous species are the chemical system's, shared by all zones).
        !! ind_sec_species holds the aqueous secondaries ONLY, starting at 1 (no primary block, unlike
        !! ind_var_act_species; see set_ind_species_react_zone). The non-aqueous blocks come from the
        !! previous-step state of the CURRENT solid chemistry (free site and adsorbed cations:
        !! solid_chemistry%conc_old, laid out as that solid chemistry) and gas chemistry
        !! (gas_chemistry%conc_old): for a displaced Lagrangian parcel, the exchanger and gas state of
        !! the cell it sits in now.
        !!
        !! The previous-step index arrays (ind_prim_species_old, ind_sec_species_old) and
        !! solid_chemistry_old, which describe the OLD zone (the one the previous-step state was
        !! equilibrated with), are NOT used here.
        !!
        !! @return conc_old_v previous-step variable activity species concentrations, size
        !!         num_var_act_species of the current reactive zone, in its column order
        !!
        function get_conc_old_var_act_new_rz(this) result(conc_old_v)
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc_old_v(:) !> previous-step variable activity concentrations, current-zone layout
            
            integer(kind=4) :: i
            
            associate(rz => this%solid_chemistry%reactive_zone, sc => this%solid_chemistry)
            allocate(conc_old_v(rz%speciation_alg%num_var_act_species))
            !> aqueous blocks: conc_old through the CURRENT index arrays (ind_sec_species starts at 1)
            do i=1,rz%speciation_alg%num_aq_prim_species
                conc_old_v(i)=this%conc_old(this%ind_prim_species(i))
            end do
            do i=1,rz%speciation_alg%num_aq_sec_var_act_species
                conc_old_v(rz%speciation_alg%num_prim_species+i)=this%conc_old(this%ind_sec_species(i))
            end do
            !> exchanger state of the previous step: the CURRENT solid chemistry's (for a displaced
            !> parcel, the cell it sits in now), laid out as that solid chemistry: minerals first, then
            !> the free site, then the adsorbed cations
            if (rz%speciation_alg%flag_cat_exch.eqv..true.) then
                conc_old_v(rz%speciation_alg%num_prim_species)=&
                    sc%conc_old(sc%mineral_zone%num_minerals+1)
                conc_old_v(rz%speciation_alg%num_prim_species+rz%speciation_alg%num_aq_sec_var_act_species+1:&
                    rz%speciation_alg%num_prim_species+rz%speciation_alg%num_aq_sec_var_act_species+&
                    rz%cat_exch_zone%num_exch_cats)=sc%conc_old(&
                    sc%mineral_zone%num_minerals+2:sc%num_solids)
            end if
            !> equilibrium variable activity gases of the previous step, from the gas chemistry
            if (associated(this%gas_chemistry)) then
                conc_old_v(rz%speciation_alg%num_aq_var_act_species+rz%cat_exch_zone%num_surf_compl+1:&
                    rz%speciation_alg%num_var_act_species-rz%gas_phase%num_gases_kin_var_act)=&
                    this%gas_chemistry%conc_old(rz%gas_phase%num_gases_eq_cst_act+1:rz%gas_phase%num_gases_eq)
            end if
            end associate
        end function
        
        !> \brief Get logarithmic activity coefficients of secondary variable activity species.
        !>
        !> This function retrieves the logarithmic activity coefficients (logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³) for all secondary
        !> variable activity species in the chemical system. Variable activity species are those whose
        !> activities change with solution composition (non-ideal behavior) and must be computed
        !> iteratively during speciation calculations.
        !>
        !> The log activity coefficients are assembled in the following order:
        !> 1. Aqueous secondary variable activity species (complexes, ions with concentration-dependent activity)
        !> 2. Equilibrium variable activity minerals (minerals with solid solution behavior)
        !> 3. Surface complexes (adsorbed species on mineral surfaces)
        !> 4. Equilibrium variable activity gases (gases with non-ideal behavior, fugacity coefficients)
        !>
        !> This is a subset of get_log_gamma2() that excludes constant activity species.
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state and activity coefficients
        !> \param[out] log_gamma2v Logarithmic activity coefficients of secondary variable activity species (allocated internally)
        !> \return log_gamma2v Vector of logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³Ã¢â€šâ€šÃ¢â€šâ„¢ÃªÅ“â‚¬ values [-] (size = num_eq_reactions)
        !>
        !> \note The size equals num_eq_reactions because each equilibrium reaction produces one
        !>       secondary variable activity species.
        !>
        !> \see get_log_gamma2, get_c2v, compute_log_act_coeffs_aq_chem
        function get_log_gamma2v(this) result(log_gamma2v)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma2v(:)
            integer(kind=4) :: i,n2v
            allocate(log_gamma2v(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            n2v=0 !> counter for secondary variable activity species
            !> Aqueous secondary variable activity species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                log_gamma2v(n2v+i)=this%log_act_coeffs(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            n2v=n2v+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            !!> Kinetic variable activity minerals
            !!> Kinetic variable activity gases
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                log_gamma2v(n2v+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            n2v=n2v+this%solid_chemistry%reactive_zone%num_minerals_var_act
            ! !> Equilibrium constant activity minerals
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                !log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=&
                !    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
                log_gamma2v(n2v+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            n2v=n2v+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma2v(n2v+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n2v=n2v+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            ! !> Ideal water
            ! !> Kinetic constant activity minerals
            
            ! !> Kinetic constant activity gases
            ! !> Equilibrium constant activity gases
        end function
        
        !> \brief Get logarithmic activity coefficients of all species in the chemical system.
        !>
        !> This function retrieves the logarithmic activity coefficients (logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³) for all species
        !> present in the reactive transport system, including aqueous, solid, and gas phases.
        !> Activity coefficients correct for non-ideal behavior in concentrated solutions and
        !> complex multi-phase systems.
        !>
        !> The log activity coefficients are assembled in the following canonical order:
        !> 1. Aqueous primary species (basis species for speciation)
        !> 2. Free sites for surface complexation (if applicable)
        !> 3. Kinetic variable activity minerals (dissolving/precipitating minerals)
        !> 4. Kinetic constant activity minerals (non-equilibrium gas species)
        !> 5. Kinetic variable activity gases
        !> 6. Kinetic constant activity gases
        !> 7. Aqueous complexes
        !> 8. Equilibrium variable activity minerals (minerals in equilibrium with solution)
        !> 9. Equilibrium constant activity minerals (log ÃŽÂ³ = 0, pure solid phases)
        !> 10. Surface complexes (adsorbed species on mineral surfaces)
        !> 11. Equilibrium variable activity gases (ideal gas behavior)
        !> 12. Equilibrium constant activity gases (ideal gas behavior)
        !> This comprehensive ordering allows reconstruction of the complete thermodynamic state
        !> for reactive transport calculations.
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state and activity coefficients
        !> \return log_gamma_vec Vector of logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³ values [-] for all species (size = num_species)
        !>
        !> \note For constant activity species (pure phases, ideal gases), log ÃŽÂ³ = 0 by definition.
        !>       Variable activity species have concentration-dependent activity coefficients.
        !>
        !> \warning This function requires that activity coefficients have been computed beforehand
        !>          using compute_log_act_coeffs_aq_chem or similar procedures.
        !>
        !> \see get_log_gamma2, get_log_gamma2v, compute_log_act_coeffs_aq_chem
        function get_log_gamma(this) result(log_gamma_vec) !> gets log_10 activity coefficients of all species
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma_vec(:)
            
            integer(kind=4) :: i, num_sp
            
            allocate(log_gamma_vec(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
            num_sp=0 !> counter for species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_gamma_vec(i)=this%log_act_coeffs(this%ind_prim_species(i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma_vec(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
                log_gamma_vec(this%solid_chemistry%reactive_zone%speciation_alg%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq-&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+i)=&
                    this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Kinetic variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
                log_gamma_vec(num_sp+i)=this%solid_chemistry%log_act_coeffs(i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
            !> Kinetic constant activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_cst_act
                log_gamma_vec(num_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_kin_cst_act
            !> Kinetic variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
                log_gamma_vec(num_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
            !> Kinetic constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_cst_act
                log_gamma_vec(num_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_cst_act
            !> Aqueous complexes
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%aq_phase%num_aq_complexes
                log_gamma_vec(num_sp+i)=this%log_act_coeffs(this%ind_sec_species(i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%chem_syst%aq_phase%num_aq_complexes
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_eq_var_act
                log_gamma_vec(num_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_eq_var_act
            !> Equilibrium constant activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
                log_gamma_vec(num_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma_vec(num_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            !> Equilibrium constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                log_gamma_vec(num_sp+i)=this%gas_chemistry%log_act_coeffs(i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        end function
        
        function get_conc(this) result(conc)
        !> gets concentrations of all species
        !! concentrations are ordered in:
            !<      primary aqueous species
            !<      free surface
            !<      kinetic variable activity minerals
            !>      kinetic variable activity gases
            !<      secondary variable activity aqueous species
            !<      equilibrium variable activity minerals
            !<      surface complexes
            !<      variable activity gases in equilibrium
            !<      ideal water
            !<      kinetic constant activity minerals
            !<      equilibrium constant activity minerals
            !<      kinetic constant activity gases
            !<      constant activity gases in equilibrium
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc(:)
            
            integer(kind=4) :: i,num_sp
            
            allocate(conc(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
            
            num_sp=0 !> counter for species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                conc(i)=this%concentrations(this%ind_var_act_species(i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                conc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=&
                    this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
                conc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)=&
                    this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Kinetic variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
                conc(num_sp+i)=this%solid_chemistry%concentrations(i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
            !> Kinetic variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
                conc(num_sp+i)=this%gas_chemistry%concentrations(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
            !> Aqueous secondary variable activity
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                conc(num_sp+i)=this%concentrations(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_eq_var_act
                conc(num_sp+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_eq_var_act+&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats !> chapuza
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                conc(num_sp+i)=this%gas_chemistry%concentrations(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            !> Ideal water
            do i=1,this%aq_phase%wat_flag
                conc(num_sp+i)=this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))
            end do
            num_sp=num_sp+this%aq_phase%wat_flag
            !> Kinetic constant activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_cst_act
                conc(num_sp+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_kin_cst_act
            !> Equilibrium constant activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
                conc(num_sp+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%num_minerals_cst_act
            !> Kinetic constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_cst_act
                conc(num_sp+i)=this%gas_chemistry%concentrations(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_cst_act
            !> Equilibrium constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                conc(num_sp+i)=this%gas_chemistry%concentrations(i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        end function
        
        !> \brief Checks for zero (or near-zero) concentrations in aqueous species
        !!
        !! \details Iterates over aqueous species concentrations and flags
        !! any that are below the squared absolute tolerance.
        !!
        !! @param zero_flag TRUE if zero concentration found
        !! @param ind_zero index of the zero concentration species
        !!
        subroutine check_zero_conc(this,zero_flag,ind_zero) !> checks zero concentrations
        implicit none
        class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry
        logical :: zero_flag !> TRUE if zero concentration, FALSE otherwise
        integer(kind=4), intent(out) :: ind_zero !> index of zero concentration
        
        integer(kind=4) :: i !> counter
        
        zero_flag=.false. !> initialize flag
        ind_zero=0 !> initialize index
        do i=1,this%aq_phase%num_species
            if (this%concentrations(i)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then
                zero_flag=.true.
                ind_zero=i
            end if
        end do
        end subroutine
        
        !> \brief Checks consistency of aqueous variable activity species concentrations
        !!
        !! \details Validates that the aqueous variable activity species concentrations
        !! are consistent with the component concentrations via the component matrix.
        !!
        !! @param conc_comp component concentrations to validate against
        !!
        subroutine check_conc_aq_var_act_species(this,conc_comp) !> checks concentration aqueous variable activity species from &
            !! components
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_comp(:)
            
            real(kind=8), allocatable :: res(:) !> residual
            
            res=conc_comp-matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,&
            this%concentrations(this%ind_var_act_species(1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)))
            if (inf_norm_vec_real(res)>this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
                error stop "Error in aqueous variable activity concentrations"
            end if
        end subroutine
        
        !> \brief Checks consistency of all variable activity species concentrations
        !!
        !! \details Validates that all variable activity species concentrations
        !! are consistent with the component concentrations via the component matrix.
        !!
        !! @param conc_nc concentration of variable activity species
        !! @param conc_comp concentration of components
        !!
        subroutine check_conc_var_act_species(this,conc_nc,conc_comp) !> checks concentration variable activity species from & 
            !! components
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_nc(:) !> concentration variable activity species
            real(kind=8), intent(in) :: conc_comp(:) !> concentration components
            
            real(kind=8), allocatable :: res(:) !> residual

            
            res=conc_comp-matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
            if (inf_norm_vec_real(res)>this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
                error stop "Error in variable activity concentrations"
            end if
        end subroutine
        
        !> \brief Checks that aqueous species activities satisfy the mass action law
        !!
        !! \details Verifies that the current activities are consistent with
        !! the equilibrium constants via the stoichiometric matrix.
        !!
        subroutine check_act_aq_species(this) !> checks activity aqueous species
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            
            real(kind=8), allocatable :: log_res(:) !> log_10(residual)
            
            !> this water's own log10 K(T) (equilibrium block), not the shared 25 C eq_cst_ref
            log_res=this%logK(1:this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,1:this%aq_phase%num_species),log10(this%activities))
            if (inf_norm_vec_real(log_res)>this%solid_chemistry%reactive_zone%CV_params%log_rel_tol) then
                print *, inf_norm_vec_real(log_res)
                error stop "Error in aqueous activities"
            end if
        end subroutine

        !> \brief Modified equilibrium constants owned by the aqueous chemistry
        !!
        !! \details Returns the reactive zone's base equilibrium constants with the
        !! constant-activity phases folded into the mass-action law:
        !!     K' = K_eq / (product of constant-activity species activities).
        !! The constant-activity phases are minerals and gases:
        !!   - Pure minerals have activity 1, so they leave K unchanged (identity).
        !!   - Constant-activity gases have activity = partial pressure P, so K' = K/P.
        !! This lives on the aqueous chemistry (not the shared reactive zone) because
        !! it needs BOTH the solid chemistry (reactive zone / mineral phases) and the
        !! gas chemistry (partial pressures), which are per-water state. Note that the
        !! constant-activity minerals' factor is 1 by definition, so only the gas
        !! partial pressures actually modify K (solid_chemistry%activities is not used
        !! here as it may hold concentration-based values during a run).
        !! The gas block position within K follows the reactive zone reaction ordering
        !! (speciation_alg%flag_comp), mirroring get_eq_csts_react_zone.
        !! \param[in] this Aqueous chemistry object (const)
        !! \return K Array of modified equilibrium constants (size num_eq_reactions)
        function get_eq_csts_aq_chem(this) result(K)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry object
            real(kind=8), allocatable :: K(:) !> equilibrium constants K_eq(T) of ALL reactions (equilibrium + kinetic)

            !> This water's equilibrium constants at its own temperature for ALL reactions:
            !> K = 10**logK over the whole logK array (equilibrium block followed by the
            !> kinetic blocks), so kinetic reactions (e.g. kinetic-mineral saturation indices)
            !> are covered too. The constant-activity phase activities (mineral a=1; gas a=P)
            !> are folded in by compute_logK_star, not here.
            K=10.0d0**this%logK
        end function

        !> \brief Builds this water's transformed (modified) equilibrium constants.
        !! \details logK_star = inv_Se_nc_2Ã‚Â·log10(K) (flag_comp=.true.), where K is this water's
        !! modified constant vector (constant-activity gas partial pressures already
        !! folded in, and already reordered by any reaction swap done in
        !! compute_speciation_alg_arrays). The structural inverse (inv_Se_2 /
        !! inv_Se_nc_2) lives on the (possibly shared) reactive zone; only these
        !! K-dependent vectors are stored per-water so that waters with different
        !! partial pressures sharing a reactive zone do not clobber each other.
        !! It reads this water's own logK directly (this%logK is already log10(K(T)) per
        !! reaction) and folds in the constant-activity phase activities in log space
        !! (mineral a=1 -> log a=0; gas a=P -> log10 K' = log10 K - log10 P) before the
        !! transform. The structural inverse (inv_Se_2 / inv_Se_nc_2) lives on the (possibly
        !! shared) reactive zone; only these K-dependent vectors are stored per-water.
        subroutine compute_logK_star(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            real(kind=8), allocatable :: logK_mod(:) !> this water's (modified) log10 equilibrium constants (equilibrium block)
            integer(kind=4) :: i,offset_gas !> loop index and gas-block offset within logK_mod
            associate(spec_alg=>this%solid_chemistry%reactive_zone%speciation_alg)
            if (spec_alg%num_eq_reactions>0) then
                !> This water's log10 equilibrium constants (equilibrium block), taken directly
                !> from this%logK (already log10(K(T)) per reaction); no 10** / log10 round-trip.
                logK_mod=this%logK(1:spec_alg%num_eq_reactions)
                !> Only logK_star is built here (inv_Se_nc_2Ã‚Â·log10(K)), from the UNMODIFIED
                !> equilibrium constants (no constant-activity gas partial-pressure fold): the gas
                !> is handled through its activity coefficient (gamma = R*T) in the mass action and
                !> its constraint P = ctot is enforced separately by the Newton residual. logK_tilde
                !> is built where it is needed (see read_wat_type_CHEPROO), not here.
                this%logK_star=matmul(spec_alg%inv_Se_nc_2,logK_mod)
            end if
            end associate
        end subroutine

        !> \brief Builds this water's transformed constant logK_tilde = inv_Se_2Ã‚Â·log10(K).
        !! \details Companion to compute_logK_star (which builds logK_star). logK_tilde is the
        !! flag_comp=.false. transform (all secondary species) consumed by the Picard speciation
        !! (compute_c2_from_c1_Picard). It is built from this water's UNMODIFIED per-reaction logK
        !! (this%logK, already log10(K(T))); the constant-activity gas partial pressure is NOT
        !! folded in (the gas is handled via its activity coefficient gamma = R*T and its
        !! constraint P = ctot is enforced by the Newton residual).
        subroutine compute_logK_tilde_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            associate(spec_alg=>this%solid_chemistry%reactive_zone%speciation_alg)
            if (spec_alg%num_eq_reactions>0) then
                this%logK_tilde=matmul(spec_alg%inv_Se_2,this%logK(1:spec_alg%num_eq_reactions))
            end if
            end associate
        end subroutine

        !> \brief Gets aqueous variable activity species concentrations after applying a permutation
        !!
        !! @param perm permutation array
        !! @return c_nc_aq permuted aqueous concentrations
        !!
        function get_conc_var_act_aq_perm(this,perm) result(c_nc_aq) !> gets concentrations of aqueous species after a permutation
           implicit none
           class(aqueous_chemistry_c), intent(in) :: this
           integer(kind=4), intent(in) :: perm(:)
           real(kind=8), allocatable :: c_nc_aq(:)
           allocate(c_nc_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
           c_nc_aq=this%concentrations(perm)
        end function
        
        !> \brief Rearranges state variables when the aqueous phase ordering changes
        !!
        !! \details Reorders concentrations, activities, log activity coefficients,
        !! and index arrays to match a new aqueous phase species ordering.
        !!
        !! @param old_aq_phase previous aqueous phase object with old species ordering
        !!
        subroutine rearrange_state_vars(this,old_aq_phase)
            implicit none
            class(aqueous_chemistry_c) :: this
            class(aq_phase_c), intent(in) :: old_aq_phase
            
            integer(kind=4) :: i,j
            real(kind=8), allocatable :: aux(:,:)
            
            if (this%aq_phase%num_species/=old_aq_phase%num_species) error stop
            !> falta comparar especies acuosas
            allocate(aux(this%aq_phase%num_species,3)) !> chapuza
            aux(:,1)=this%concentrations
            aux(:,2)=this%activities
            aux(:,3)=this%log_act_coeffs
            i=1
            j=1
            do 
                if (this%aq_phase%aq_species(i)%name==old_aq_phase%aq_species(j)%name) then
                    this%concentrations(i)=aux(j,1)
                    this%activities(i)=aux(j,2)                                   !< Restore activity value from temporary storage
                    this%log_act_coeffs(i)=aux(j,3)                        !< Restore log activity coefficient from temporary
                    if (i<this%aq_phase%num_species) then                  !< If not at last species
                        i=i+1                                               !< Move to next species in current phase
                        j=1                                                 !< Reset old phase index to beginning
                    else                                                    !< If at last species
                        exit                                                !< Exit the main loop
                    end if
                else if (j<old_aq_phase%num_species) then                   !< If more species in old phase to check
                    j=j+1                                                   !< Move to next species in old phase
                else if (i<this%aq_phase%num_species) then                  !< If more species in current phase
                    i=i+1                                                   !< Move to next species in current phase
                    j=1                                                     !< Reset old phase index to beginning
                else                                                        !< If no more species to process
                    exit                                                    !< Exit the main loop
                end if
            end do
        end subroutine
        
        !> Validate Jacobian computation for secondary variable activity species (aqueous version)
        !! Mathematical Context: Checks consistency of analytical Jacobian dc2v/dc1 with activity coefficient derivatives
        !! This verifies the chain rule application in speciation calculations with activity coefficients
        subroutine check_dc2v_dc1_aq(this,c2v,dc2v_dc1,log_Jacobian_act_coeffs)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), intent(in) :: c2v(:)                            !< Secondary variable activity concentrations
            real(kind=8), intent(in) :: dc2v_dc1(:,:)                     !< Jacobian of secondary w.r.t primary concentrations
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)       !< Jacobian of log activity coefficients
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)  !< Matrix computation variables
            type(diag_matrix_c) :: c1_diag,c2v_inv_diag                   !< Diagonal matrices for scaling
            
            !> Create diagonal matrices for numerical scaling
            call c1_diag%set_diag_matrix(this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            call c2v_inv_diag%set_diag_matrix(1d0/c2v)                   !< Inverse concentration diagonal matrix
            !> Compute left-hand side: activity coefficient contribution
            lhs_1=c2v_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))+log_Jacobian_act_coeffs(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
                                                                            !< Complex matrix expression for activity effects
            !> Multiply by concentration scaling
            lhs_2=c1_diag%prod_mat_diag_mat(dc2v_dc1)                      !< Scale Jacobian by primary concentrations
            lhs=matmul(lhs_1,lhs_2)                                         !< Complete left-hand side computation
            !> Compute right-hand side: expected value from stoichiometry
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(&
            1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+&
            id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
                                                                            !< Right-hand side from theoretical derivation
            !> Check consistency within tolerance
            res=lhs-rhs                                                     !< Compute residual between analytical and expected
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_rel_tol) then
                print *, norm_mat_inf(res)                                  !< Print error magnitude for debugging
                error stop "Error in dc2v_dc1"                            !< Stop execution if validation fails
            end if
        end subroutine
        
        !> Validate Jacobian computation for secondary variable activity species (general version)
        !! Mathematical Context: Similar to aqueous version but uses external concentration arrays
        !! This allows validation of Jacobians computed with different concentration sets
        subroutine check_dc2v_dc1(this,c1,c2v,dc2v_dc1,log_Jacobian_act_coeffs)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations (external array - "chapuza")
            real(kind=8), intent(in) :: c2v(:)                            !< Secondary variable activity concentrations ("chapuza")
            real(kind=8), intent(in) :: dc2v_dc1(:,:)                     !< Jacobian matrix to validate
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)       !< Log activity coefficient Jacobian
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)  !< Matrix computation variables
            type(diag_matrix_c) :: c1_diag,c2_inv_diag                     !< Diagonal matrices for scaling
            
            !> Create diagonal matrices using external concentration arrays
            call c1_diag%set_diag_matrix(c1)                               !< Primary concentration diagonal matrix
            call c2_inv_diag%set_diag_matrix(1d0/c2v)                     !< Inverse secondary concentration diagonal
            !> Compute left-hand side with activity coefficient effects
            lhs_1=c2_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))+&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
                                                                            !< Activity coefficient matrix transformation
            lhs_2=c1_diag%prod_mat_diag_mat(dc2v_dc1)                      !< Scale by primary concentrations
            lhs=matmul(lhs_1,lhs_2)                                         !< Complete left-hand side
            !> Compute theoretical right-hand side
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+id_matrix(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-log_Jacobian_act_coeffs(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
                                                                            !< Expected value from analytical derivation
            !> Validate within numerical tolerance
            res=lhs-rhs                                                     !< Compute validation residual
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_rel_tol) then
                print *, norm_mat_inf(res)                                  !< Print error for debugging
                error stop "Error in dc2v_dc1"                            !< Terminate if validation fails
            end if
        end subroutine
        
        !> Validate Jacobian computation for all secondary species (including constant activity)
        !! Mathematical Context: Extends validation to include secondary species with constant activities
        !! This is the most general validation covering the complete speciation Jacobian
        subroutine check_dc2_dc1(this,c1,c2,dc2_dc1,log_Jacobian_act_coeffs)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c), intent(in) :: this                 !< Aqueous chemistry object (input only)
            real(kind=8), intent(in) :: c1(:)                              !< Primary concentrations (external - "chapuza")
            real(kind=8), intent(in) :: c2(:)                              !< All secondary concentrations ("chapuza")
            real(kind=8), intent(in) :: dc2_dc1(:,:)                       !< Complete secondary Jacobian matrix
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)       !< Log activity coefficient Jacobian
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)
            type(diag_matrix_c) :: c1_diag,c2_inv_diag
            
            call c1_diag%set_diag_matrix(c1(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            call c2_inv_diag%set_diag_matrix(1d0/c2)
            lhs_1=c2_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species))+&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species))
            lhs_2=c1_diag%prod_mat_diag_mat(dc2_dc1)
            lhs=matmul(lhs_1,lhs_2)
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_Jacobian_act_coeffs(&
            1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+&
            id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            res=lhs-rhs
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_rel_tol) then
                print *, norm_mat_inf(res)
                error stop "Error in dc2_dc1"
            end if
        end subroutine
        
        !> \brief Computes the logarithmic activity coefficient of water
        !!
        !! \details Calculates log(gamma_w) from the ratio of the activity
        !! to the concentration of water: log(gamma) = log10(a) - log10(c).
        !!
        subroutine compute_log_act_coeff_wat(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (this%aq_phase%ind_wat>0) then
                this%log_act_coeffs(this%indices_aq_species(this%aq_phase%ind_wat))=&
                    log10(this%activities(this%indices_aq_species(this%aq_phase%ind_wat)))-&
                    log10(this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat)))
            end if
        end subroutine
        
        !> \brief Allocates the per-water saturation indices: one entry per mineral of the chemical
        !!        system (kinetic first, then equilibrium, as chem_syst%minerals), initialised to 1
        !!        (equilibrium: Omega = 1, no kinetic driving force).
        !! \details Idempotent: an array of the right size is kept (its values untouched); a stale
        !!          one of another size (chemical system changed) is reallocated and reset to 1.
        !!          Requires the solid chemistry, hence the chemical system, to be associated.
        subroutine allocate_saturation_mins(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: n_min !> number of minerals of the chemical system (kinetic + equilibrium)
            if (.not. associated(this%solid_chemistry)) error stop "allocate_saturation_mins: solid chemistry not associated"
            !> sized by the CHEMICAL SYSTEM, not by the water's zone: the array is indexed by the mineral's
            !> position in chem_syst%minerals (compute_saturation_mins loops over all of them)
            n_min=this%solid_chemistry%reactive_zone%chem_syst%num_minerals
            if (allocated(this%saturation_mins)) then
                if (size(this%saturation_mins)/=n_min) deallocate(this%saturation_mins)
            end if
            if (.not. allocated(this%saturation_mins)) then
                allocate(this%saturation_mins(n_min))
                this%saturation_mins=1d0 !> never evaluated = equilibrium (Omega=1, no kinetic driving force)
            end if
        end subroutine allocate_saturation_mins

        function compute_saturation_kin_min(this,ind_min) result(saturation)
        !> This function is meant to be used for minerals NOT in equilibrium.
        !> Besides returning Omega = IAP/K(T), it stores it in this%saturation_mins(ind_min), the per-water
        !> record of the saturation indices of all minerals (kinetic minerals come first in
        !> chem_syst%minerals, so the kinetic index ind_min is also the mineral's position there).
        !> [SAT-DBG] When dbg_kin_mineral (kin_mineral_m) is armed by the solver — first target water, first
        !> dbg_kin_mineral_max_step time step(s) of the simulation (default: the first step only) — every factor
        !> of Omega (species, their concentration, log gamma, activity and a^(-nu), the IAP, the logK(T) used and
        !> its index) is written to the <root>_kin_mineral_dbg.out file (standard output if not open), one block
        !> per call.
            use kin_mineral_m, only: dbg_kin_mineral, dbg_kin_mineral_unit, dbg_kin_mineral_step !> [SAT-DBG] trace switch, file unit and current step
            use iso_fortran_env, only: output_unit !> fallback of the [SAT-DBG] trace when its file is not open
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: ind_min !> index of kinetic mineral in chemical system
            real(kind=8) :: saturation !> saturation index
            
            real(kind=8) :: IAP
            integer(kind=4) :: i,ind_kin_react
            integer(kind=4) :: ind_logK !> position of the mineral's logK(T) in this%logK (equilibrium, linear and redox kinetic reactions first)
            integer(kind=4) :: ind_aq_ph !> [SAT-DBG] aqueous phase index of species i of the reaction
            integer(kind=4) :: ind_st !> [SAT-DBG] index of species i in this water's state arrays (concentrations/activities)
            integer(kind=4) :: iu !> [SAT-DBG] output unit of the trace
            integer(kind=4), save :: ncall_dbg=0 !> [SAT-DBG] number of traced calls so far
            real(kind=8) :: nu_i !> [SAT-DBG] stoichiometric coefficient of species i
            real(kind=8) :: lg_i !> [SAT-DBG] log10 activity coefficient of species i (0 if not stored)
            logical :: trace !> [SAT-DBG] write this call?

            call this%allocate_saturation_mins()
            associate(react => this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(ind_min))
            ind_logK=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions+&
                this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts+ind_min
            trace=dbg_kin_mineral !> armed by the solver for the first target water during the first time step(s) only
            if (trace) then
                ncall_dbg=ncall_dbg+1
                iu=output_unit
                if (dbg_kin_mineral_unit/=-1) iu=dbg_kin_mineral_unit
                write(iu,'(A)') '[SAT-DBG] ======================================================================'
                write(iu,'(A,A,A,I0,A,I0,A,I0)') '[SAT-DBG] compute_saturation_kin_min: mineral ', trim(react%name), &
                    ' (ind_min = ', ind_min, '), time step ', dbg_kin_mineral_step, ', traced call ', ncall_dbg
                write(iu,'(A,A,A,I0,A,ES14.6,A)') '[SAT-DBG]   water "', trim(this%name), '" (id ', this%id, '), temp = ', &
                    this%temp, ' K'
                write(iu,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '[SAT-DBG]   num_aq_rk = ', react%num_aq_rk, &
                    ' ; logK index = num_eq_reactions (', this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions, &
                    ') + num_lin_kin (', this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts, &
                    ') + num_redox_kin (', this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts, &
                    ') + ind_min (', ind_min, ') = ', ind_logK
                write(iu,'(A)') '[SAT-DBG]     i  aq_ph  state  species              nu              conc       log_gamma'//&
                    '        activity         a^(-nu)'
            end if
            IAP=1d0
            do i=1,react%num_aq_rk
                ind_aq_ph=react%indices_aq_phase(i)
                ind_st=this%indices_aq_species(ind_aq_ph)
                nu_i=react%stoichiometry(react%indices_react_species(i))
                IAP=IAP*this%activities(ind_st)**(-nu_i)
                if (trace) then
                    lg_i=0d0
                    if (allocated(this%log_act_coeffs)) lg_i=this%log_act_coeffs(ind_st)
                    write(iu,'(A,I5,2X,I5,2X,I5,2X,A20,1X,F8.3,4(1X,ES15.6))') '[SAT-DBG] ', i, ind_aq_ph, ind_st, &
                        adjustl(this%aq_phase%aq_species(ind_aq_ph)%name), nu_i, this%concentrations(ind_st), lg_i, &
                        this%activities(ind_st), this%activities(ind_st)**(-nu_i)
                end if
            end do
            saturation=IAP*(1d1**this%logK(ind_logK))
            this%saturation_mins(ind_min)=saturation
            if (trace) then
                write(iu,'(A,ES15.6)') '[SAT-DBG]   IAP = prod_i a_i^(-nu_i)          = ', IAP
                write(iu,'(A,ES15.6,A,ES15.6)') '[SAT-DBG]   logK(T) = ', this%logK(ind_logK), ', 10**logK(T) = ', &
                    1d1**this%logK(ind_logK)
                write(iu,'(A,ES15.6,A,I0,A)') '[SAT-DBG]   Omega = IAP * 10**logK(T)           = ', saturation, &
                    '  (stored in saturation_mins(', ind_min, '))'
                write(iu,'(A)') '[SAT-DBG] ======================================================================'
            end if
            end associate
        end function

        !> \brief Saturation index Omega = IAP/K(T) of an EQUILIBRIUM mineral of the chemical system,
        !!        stored in this%saturation_mins(ind_min) as compute_saturation_kin_min does for kinetic ones.
        !! \details Same convention as the kinetic minerals (both reaction sets are read alike from the
        !!          database: mineral stoichiometry +1, aqueous species -nu, K = 10^-logK):
        !!            Omega = 10^logK(T) * prod_i a_i^(-nu_i)
        !!          over the aqueous species of the mineral's equilibrium reaction, with the water's own
        !!          activities and temperature. For a mineral of the water's reactive zone this is 1 up
        !!          to the speciation tolerance; for an absent mineral it tells whether the water would
        !!          dissolve (<1) or precipitate (>1) it. The reaction is located in chem_syst%eq_reacts
        !!          by its product (the mineral, last species); its aqueous species are mapped to the
        !!          aqueous phase by name.
        function compute_saturation_eq_min(this,ind_min) result(saturation)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4), intent(in) :: ind_min !> position of the (equilibrium) mineral in chem_syst%minerals (> num_minerals_kin)
            real(kind=8) :: saturation !> saturation index Omega [-]
            
            real(kind=8) :: IAP !> ion activity product prod_i a_i^(-nu_i)
            integer(kind=4) :: i,j,ind_react,sp_ind,aq_sp_ind
            logical :: flag

            call this%allocate_saturation_mins()
            associate(cs => this%solid_chemistry%reactive_zone%chem_syst)
            !> equilibrium reaction whose product (last species) is this mineral
            ind_react=0
            do j=1,cs%speciation_alg%num_eq_reactions
                if (cs%eq_reacts(j)%react_type/=2) cycle
                sp_ind=cs%eq_reacts(j)%species_ind(cs%eq_reacts(j)%num_species)
                if (cs%species(sp_ind)%name==cs%minerals(ind_min)%name) then
                    ind_react=j
                    exit
                end if
            end do
            if (ind_react==0) then
                print *, "ERROR compute_saturation_eq_min: no equilibrium reaction with product ", trim(cs%minerals(ind_min)%name)
                error stop "compute_saturation_eq_min: equilibrium mineral without equilibrium reaction in the chemical system"
            end if
            IAP=1d0
            do i=1,cs%eq_reacts(ind_react)%num_species-1 !> reactants: the aqueous species (the mineral is the last species)
                call this%aq_phase%is_species_in_aq_phase(cs%species(cs%eq_reacts(ind_react)%species_ind(i)),flag,aq_sp_ind)
                if (.not. flag) then
                    print *, "ERROR compute_saturation_eq_min: species ", &
                        trim(cs%species(cs%eq_reacts(ind_react)%species_ind(i))%name), " of the reaction of mineral ", &
                        trim(cs%minerals(ind_min)%name), " is not in the aqueous phase"
                    error stop "compute_saturation_eq_min: non-aqueous reactant in a mineral equilibrium reaction"
                end if
                IAP=IAP*this%activities(this%indices_aq_species(aq_sp_ind))**(-cs%eq_reacts(ind_react)%stoichiometry(i))
            end do
            saturation=IAP*(1d1**logK_react_T(cs%eq_reacts(ind_react),this%temp))
            end associate
            this%saturation_mins(ind_min)=saturation
        end function compute_saturation_eq_min

        !> \brief Saturation indices of ALL minerals of the chemical system (kinetic, then equilibrium)
        !!        at the water's current state, stored in this%saturation_mins. Called at initialisation and
        !!        at the end of each time step (with the other derived quantities, pH etc.), so the
        !!        output reports the end-of-step state for every mineral, present in the water's
        !!        reactive zone or not. Requires the per-water logK (compute_logK_PHREEQC_aq_chem).
        subroutine compute_saturation_mins(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: i
            real(kind=8) :: sat
            if (.not. allocated(this%logK) .or. .not. allocated(this%activities)) return !> water not speciated yet
            associate(cs => this%solid_chemistry%reactive_zone%chem_syst)
            do i=1,cs%num_minerals_kin
                sat=this%compute_saturation_kin_min(i)
            end do
            do i=cs%num_minerals_kin+1,cs%num_minerals
                sat=this%compute_saturation_eq_min(i)
            end do
            end associate
        end subroutine compute_saturation_mins
        
        !> \brief Computes equilibrium gas partial pressures from the mass-action law.
        !!
        !! \details For each equilibrium gas the partial pressure (its activity) is
        !! obtained directly from the equilibrium mass-action law using the already
        !! speciated aqueous activities and this water's logK, instead of the ideal-gas
        !! conversion P = c*R*T/V (which mishandles the gas fugacity coefficient and
        !! yields partial pressures that are ~2 orders of magnitude too high).
        !! For the equilibrium reaction r that defines gas g, the mass-action law reads
        !!     sum_i S_ri * log10(a_i) = log10(K_r)     (sum over ALL species i)
        !! Only aqueous species and the gas itself appear in the usual gas reactions
        !! (e.g. co2(g) + h2o = h+ + hco3-), so isolating the gas term (coefficient s_g):
        !!     log10(P_g) = ( log10(K_r) - sum_aq S_rj*log10(a_j) ) / s_g
        !! The gas column in the reactive-zone stoichiometric matrix is
        !! ind_gases_stoich_mat(g); a gas is a secondary species that appears in exactly
        !! one equilibrium reaction, so its reaction row is the unique nonzero entry of
        !! that column. The aqueous columns 1:num_aq_species align with this%activities
        !! (same convention as check_act_aq_species). Results overwrite
        !! gas_chemistry%activities (the reported partial pressures).
        subroutine compute_partial_pressures_eq(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: g,r,r_g,col_g,j,n_aq !> loop/index variables
            real(kind=8) :: aq_contrib,s_g !> aqueous mass-action contribution and gas stoich coefficient
            real(kind=8), parameter :: tol=1d-12 !> tolerance to treat a stoichiometric coefficient as nonzero
            if (.not. associated(this%gas_chemistry)) return
            associate(rz=>this%solid_chemistry%reactive_zone)
                if (rz%gas_phase%num_gases_eq<1) return
                n_aq=this%aq_phase%num_species
                do g=1,rz%gas_phase%num_gases_eq
                    col_g=rz%ind_gases_stoich_mat(g) !> column of gas g in the stoichiometric matrix
                    !> Locate the (unique) equilibrium reaction that defines this gas
                    r_g=0
                    do r=1,rz%speciation_alg%num_eq_reactions
                        if (abs(rz%stoich_mat(r,col_g))>tol) then
                            r_g=r
                            exit
                        end if
                    end do
                    if (r_g==0) cycle !> gas not found in any equilibrium reaction (should not happen)
                    s_g=rz%stoich_mat(r_g,col_g) !> gas stoichiometric coefficient in its own reaction
                    !> Aqueous contribution sum_aq S_rj*log10(a_j); skip zero coefficients so
                    !> species with vanishing activity (log10 -> -Inf) never enter the sum.
                    aq_contrib=0d0
                    do j=1,n_aq
                        if (abs(rz%stoich_mat(r_g,j))>tol) then
                            aq_contrib=aq_contrib+rz%stoich_mat(r_g,j)*log10(this%activities(j))
                        end if
                    end do
                    !> log10(P_g) = ( log10(K_r) - aqueous contribution ) / s_g
                    this%gas_chemistry%activities(g)=1d1**((this%logK(r_g)-aq_contrib)/s_g)
                end do
            end associate
        end subroutine
        
     
        !> \brief Computes concentrations from component concentration for a binary system
        !!
        !! \details Uses the analytical solution of the quadratic equation
        !! arising from the mass action law in a binary system.
        !!
        !! @param u component concentration
        !! @param K equilibrium constant
        !!
        subroutine compute_c_from_u_bin_syst(this,u,K)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: u
            real(kind=8), intent(in) :: K
            this%concentrations(1)=(u+sqrt(u**2+4*K))/2d0
            this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1)=(-u+sqrt(u**2+4*K))/2d0
        end subroutine
        
        !> \brief Initialises primary species concentrations for iterative methods
        !!
        !! \details Uses linear extrapolation from two previous time steps:
        !! c = (1+param)*c_old - param*c_old_old. Falls back to c_old if
        !! extrapolation yields a negative value.
        !!
        !! @param conc_old_old concentrations from two time steps ago
        !! @param conc_old concentrations from previous time step
        !! @param param extrapolation parameter in [0,1]
        !!
        subroutine initialise_conc_nc_iterative_method_bis(this,conc_old_old,conc_old,param)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_old_old(:)
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: param
            
            integer(kind=4) :: i
            
            if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                if ((1d0+param)*conc_old(i)-param*conc_old_old(i)<0d0) then
                    this%concentrations(this%ind_var_act_species(i))=conc_old(i)
                else
                    this%concentrations(this%ind_var_act_species(i))=(1d0+param)*conc_old(i)-param*conc_old_old(i)
                end if
            end do
        end subroutine
        
     
        !> \brief Sets activities of aqueous species
        !!
        !! \details If no argument provided, sets ideal activities:
        !! dissolved solids activity = concentration, water activity = 1.
        !!
        !! @param act optional activities array
        !!
        subroutine set_act_aq_species(this,act)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: act(:)
            integer(kind=4) :: i
            if (present(act)) then
                do i=1,this%aq_phase%num_species
                    this%activities(i)=act(i)
                end do
            else !> ideal activities by default
                this%activities(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                    this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))
                this%activities(this%indices_aq_species(this%aq_phase%ind_wat))=1d0 !> activity of ideal water
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=1d0/18d-3 !> concentration of ideal water in mol/m3
            end if
        end subroutine
        
        !> Allocate memory for indices of aqueous species
        subroutine allocate_indices_aq_species_aq_chem(this)
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%indices_aq_species)) deallocate(this%indices_aq_species) !< Deallocate if already allocated
            allocate(this%indices_aq_species(this%aq_phase%num_species))   !< Allocate array for species indices
        end subroutine
        
        !> Allocate memory for indices of aqueous phase mapping
        subroutine allocate_indices_aq_phase_aq_chem(this)
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (allocated(this%indices_aq_phase)) deallocate(this%indices_aq_phase) !< Deallocate if already allocated
            allocate(this%indices_aq_phase(this%aq_phase%num_species))     !< Allocate array for phase indices
        end subroutine
        
        !> Set default indices for aqueous species and phase mapping
        subroutine set_indices_aq_species_aq_chem(this)
            
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            
            integer(kind=4) :: i                                           !< Loop index for species iteration
            
            call this%allocate_indices_aq_species_aq_chem()                !< Allocate memory for species indices
            call this%allocate_indices_aq_phase_aq_chem()                  !< Allocate memory for phase indices
            do i=1,this%aq_phase%num_species                               !< Loop over all aqueous species
                this%indices_aq_species(i)=i                               !< Set default species index (1:1 mapping)
                this%indices_aq_phase(i)=i                                 !< Set default phase index (1:1 mapping)
            end do
        end subroutine
    

        !> Compute Jacobian of ionic strength with respect to primary concentrations
        !! Mathematical Context: I = 1/2 * ÃŽÂ£(ci * ziÃ‚Â²) where I is ionic strength, ci is concentration, zi is charge
        !! The derivative dI/dci = 1/2 * ziÃ‚Â² + ÃŽÂ£(dcsec/dci * zsecÃ‚Â²) for secondary species dependencies
        subroutine compute_dI_dc1(this,dc2aq_dc1,dI_dc1)
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: dc2aq_dc1(:,:)                     !< Jacobian of secondary aqueous concentrations w.r.t primary
            real(kind=8), intent(out) :: dI_dc1(:)                         !< Jacobian of ionic strength w.r.t primary concentrations
                                                                            !! (must be already allocated)
            
            integer(kind=4) :: i                                           !< Loop index for primary species
            
            !> For each primary species, compute dI/dci = 1/2*ziÃ‚Â² + ÃŽÂ£(dcsec/dci * zsecÃ‚Â²)
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                dI_dc1(i)=5d-1*this%aq_phase%aq_species(i)%valence**2 + &  !< Direct contribution: 1/2 * ziÃ‚Â²
                dot_product(this%solid_chemistry%reactive_zone%chem_syst%z2(&
                this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species),dc2aq_dc1(:,i))
                                                                            !< Secondary species contribution via chain rule
            end do

        end subroutine

       
        !> Compute Jacobian of Newton residual function for reactive mixing with Euler Fully Implicit (standard approach)
        !! Mathematical Context: Similar to lumped version but includes explicit mixing ratio for reaction rates
        !! This allows for different time stepping strategies in transport vs. chemistry
        subroutine compute_dfk_dc_aq_EfI(this,drk_dc,Delta_t,theta,mix_ratio_Rk,dfk_dc)
            !> Arguments
                class(aqueous_chemistry_c), intent(in) :: this             !< Aqueous chemistry object (polymorphic)
                real(kind=8), intent(in) :: drk_dc(:,:)                    !< Jacobian of kinetic reaction rates w.r.t concentrations
                !real(kind=8), intent(in) :: porosity                      !< Porosity (commented - in solid chemistry)
                real(kind=8), intent(in) :: Delta_t                        !< Time step size for temporal discretization
                real(kind=8), intent(in) :: theta                          !< Time weighting factor (0=explicit, 1=implicit)
                real(kind=8), intent(in) :: mix_ratio_Rk                   !< Mixing ratio for kinetic reaction rates
                real(kind=8), intent(out) :: dfk_dc(:,:)                   !< Jacobian Newton residual w.r.t concentrations (pre-allocated)
            !> Process
                !> Compute Newton residual Jacobian with mixing ratio: df/dc = I - ÃŽÂ¸*ÃŽâ€t*mix_ratio*U*S*dR/dc
                dfk_dc=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)-theta*Delta_t*mix_ratio_Rk*&
                    matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc)  !< Apply stoichiometric transformation with mixing
            end subroutine

        !> Compute component concentrations after mixing (transport step)
        !! Mathematical Context: u_mix represents component concentrations after advection-dispersion
        !! but before chemical reactions. This is the coupling point between transport and chemistry.
        !! Formula: u_mix = U * c_mix where U is the component matrix
        function compute_u_mix(this,c_mix) result(u_mix)
            implicit none                                                           !< No implicit variable declarations
        !> Arguments
            class(aqueous_chemistry_c), intent(in) :: this                         !< Aqueous chemistry object (input only)
            real(kind=8), intent(in) :: c_mix(:)                                 !< Concentration of "mobile" species after mixing
            real(kind=8), allocatable :: u_mix(:)                                !< Component concentrations after mixing (output)
        !> Variables
            integer(kind=4) :: i                                                   !< Index for aqueous components (unused)
            integer(kind=4) :: j                                                   !< Index for mixing waters (unused)
        !> Process
            !> Transform species concentrations to component concentrations using component matrix
            u_mix=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_mix)
                                                                                    !< Matrix-vector product: U * c_mix = u_mix
        end function

        !> Write a legend for the [CELL-DBG] output.
        !! \details Prints (1) the component definitions as linear combinations of the
        !! variable-activity species, in the same style as u_init.dat, and (2) the reaction
        !! names in the same order as the reaction-rate value rows (aqueous equilibrium,
        !! mineral+exchange equilibrium, aqueous kinetic, mineral kinetic). Names are listed
        !! in the current speciation/array order (values are NOT reordered).
        subroutine write_celldbg_legend(this,iunit)
            implicit none                                                           !< No implicit variable declarations
        !> Arguments
            class(aqueous_chemistry_c), intent(in) :: this                         !< Aqueous chemistry object (input only)
            integer(kind=4), intent(in) :: iunit                                   !< Output unit of the CELL-DBG file
        !> Variables
            integer(kind=4) :: ci                                                  !< Component (comp_mat row) index
            integer(kind=4) :: j                                                   !< Species / reaction index
            integer(kind=4) :: n_c                                                 !< Number of components
            integer(kind=4) :: n_v                                                 !< Number of variable-activity species
            integer(kind=4) :: n_aq_prim                                           !< Number of aqueous primary species (secondary aqueous species follow)
            integer(kind=4) :: n_min                                               !< Number of equilibrium minerals
            integer(kind=4) :: n_exch                                              !< Number of exchangeable cations
            integer(kind=4) :: n_aqeq                                              !< Number of aqueous equilibrium reactions
            integer(kind=4) :: n_lin                                               !< Number of linear kinetic reactions
            integer(kind=4) :: n_redox                                             !< Number of redox kinetic reactions
            integer(kind=4) :: n_minkin                                            !< Number of kinetic minerals
            real(kind=8) :: coef                                                   !< Component-matrix coefficient
            real(kind=8), parameter :: tol_coef=1d-12                              !< Tolerance to treat a coefficient as zero
            logical :: first                                                       !< Flag for the first term of a combination
            character(len=1024) :: expr                                            !< Accumulated linear-combination string
            character(len=64) :: term                                             !< Single term of the linear combination
            character(len=64) :: coef_str                                         !< Formatted coefficient
            character(len=:), allocatable :: sp_name                              !< Variable-activity species name
        !> Process
            associate(sc => this%solid_chemistry, rz => this%solid_chemistry%reactive_zone, &
                      cs => this%solid_chemistry%reactive_zone%chem_syst)
            n_c=rz%speciation_alg%num_prim_species                                 !< comp_mat rows (components)
            n_v=rz%speciation_alg%num_var_act_species                              !< comp_mat columns (ALL variable-activity species: aqueous + exchanger)
            n_aq_prim=rz%speciation_alg%num_aq_prim_species                        !< secondary aqueous species occupy var-act indices n_aq_prim+1..n_v
            !> Component definitions as linear combinations of variable-activity species (u_init.dat style)
            write(iunit,'(A)') "[CELL-DBG]   components (u_i = combination of variable-activity species):"
            do ci=1,n_c
                expr=''
                first=.true.
                do j=1,n_v
                    coef=rz%speciation_alg%comp_mat(ci,j)
                    if (abs(coef)<tol_coef) cycle
                    !> comp_mat column j maps to chem_syst species rz%ind_var_act_species(j);
                    !> this spans ALL variable-activity species (aqueous, free exchange site
                    !> and surface complexes), so the exchanger component prints its x-/surface
                    !> terms instead of collapsing to 0.
                    sp_name=trim(cs%species(rz%ind_var_act_species(j))%name)
                    if (first) then
                        if (coef>0d0) then
                            if (abs(coef-1d0)<tol_coef) then
                                term=sp_name
                            else
                                write(coef_str,'(F0.1)') coef
                                term=trim(adjustl(coef_str))//' '//sp_name
                            end if
                        else
                            if (abs(coef+1d0)<tol_coef) then
                                term='- '//sp_name
                            else
                                write(coef_str,'(F0.1)') abs(coef)
                                term='- '//trim(adjustl(coef_str))//' '//sp_name
                            end if
                        end if
                        first=.false.
                    else
                        if (coef>0d0) then
                            if (abs(coef-1d0)<tol_coef) then
                                term=' + '//sp_name
                            else
                                write(coef_str,'(F0.1)') coef
                                term=' + '//trim(adjustl(coef_str))//' '//sp_name
                            end if
                        else
                            if (abs(coef+1d0)<tol_coef) then
                                term=' - '//sp_name
                            else
                                write(coef_str,'(F0.1)') abs(coef)
                                term=' - '//trim(adjustl(coef_str))//' '//sp_name
                            end if
                        end if
                    end if
                    expr=trim(expr)//trim(term)
                end do
                if (first) expr='0'
                write(iunit,'(A,I0,A,A)') "[CELL-DBG]     u", ci, " = ", trim(expr)
            end do
            !> Reaction names, one per line, in the same order as the reaction-rate value rows below
            !> Each aqueous equilibrium reaction eliminates the secondary aqueous species at var-act index n_aq_prim+j
            n_aqeq=cs%num_aq_eq_reacts
            if (n_aqeq>0) then
                write(iunit,'(A)') "[CELL-DBG]   aq eq reactions (secondary species eliminated, aq re_mean order):"
                do j=1,n_aqeq
                    if (n_aq_prim+j<=n_v) then
                        write(iunit,'(A,I0,A,A)') "[CELL-DBG]     re", j, " = ", &
                            trim(this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(n_aq_prim+j)))%name)
                    else
                        write(iunit,'(A,I0,A,A)') "[CELL-DBG]     re", j, " = ", trim(cs%eq_reacts(j)%name)
                    end if
                end do
            end if
            n_min=rz%num_minerals
            n_exch=rz%cat_exch_zone%num_exch_cats
            !> Mineral equilibrium reactions (first block of the min/exch re_mean & amnt_eq value rows)
            if (n_min>0) then
                write(iunit,'(A)') "[CELL-DBG]   mineral eq reactions (min re_mean / min amnt_eq order, entries 1..n_min):"
                do j=1,n_min
                    write(iunit,'(A,I0,A,A)') "[CELL-DBG]     me", j, " = ", &
                        trim(cs%minerals(rz%ind_mins_chem_syst(j))%name)
                end do
            end if
            !> Exchange (cation-exchange) equilibrium reactions (they follow the minerals in the
            !> same re_mean & amnt_eq value rows, at entries n_min+1..n_min+n_exch). surf_compl(1) is
            !> the free site (exchanger master/primary species, not a reaction); the n_exch exchange
            !> reactions produce the adsorbed cations, which are surf_compl(2:1+n_exch).
            if (n_exch>0) then
                write(iunit,'(A)') "[CELL-DBG]   exchange eq reactions:"
                do j=1,n_exch
                    write(iunit,'(A,I0,A,A)') "[CELL-DBG]     me", n_min+j, " = ", &
                        trim(rz%cat_exch_zone%surf_compl(1+j)%name)
                end do
            end if
            n_lin=cs%num_lin_kin_reacts
            n_redox=cs%num_redox_kin_reacts
            if (n_lin+n_redox>0) then
                write(iunit,'(A)') "[CELL-DBG]   aq kin reactions (aq rk_mean/new/old order):"
                do j=1,n_lin
                    write(iunit,'(A,I0,A,A)') "[CELL-DBG]     rk", j, " = ", &
                        trim(cs%lin_kin_reacts(j)%name)
                end do
                do j=1,n_redox
                    write(iunit,'(A,I0,A,A)') "[CELL-DBG]     rk", n_lin+j, " = ", &
                        trim(cs%redox_kin_reacts(j)%name)
                end do
            end if
            n_minkin=sc%mineral_zone%num_minerals_kin
            if (n_minkin>0) then
                write(iunit,'(A)') "[CELL-DBG]   min kin reactions (min rk_mean/new/old order):"
                do j=1,n_minkin
                    write(iunit,'(A,I0,A,A)') "[CELL-DBG]     mrk", j, " = ", &
                        trim(cs%minerals(sc%mineral_zone%ind_min_chem_syst(j))%name)
                end do
            end if
            end associate
        end subroutine

        !> Write one cell's [CELL-DBG] trace block (step header, species, legend on k==1,
        !! c/u mixed & conservative & new values, reaction contributions and rates).
        !! Centralises the per-cell debug output so every solver emits identical rows.
        !! c_mix is the transported mix (lumped solvers pass their c_mix/c_hat). c_tilde
        !! (conservative part) and c_hat (reacted mix, non-lumped solvers) are optional and,
        !! when present, add the c_tilde/u_tilde and c_hat/u_hat rows respectively.
        subroutine write_celldbg_cell(this,iunit,k,twi,n_v,c_mix,c_tilde,c_hat)
            implicit none                                                           !< No implicit variable declarations
        !> Arguments
            class(aqueous_chemistry_c), intent(in) :: this                         !< Aqueous chemistry object (input only)
            integer(kind=4), intent(in) :: iunit                                   !< Output unit of the CELL-DBG file
            integer(kind=4), intent(in) :: k                                       !< Current time step index
            integer(kind=4), intent(in) :: twi                                    !< Target water index
            integer(kind=4), intent(in) :: n_v                                     !< Number of variable-activity species
            real(kind=8), intent(in) :: c_mix(:)                                   !< Transported (mixed) var-act species conc
            real(kind=8), intent(in) :: c_tilde(:)                                 !< Conservative part of the mix
            real(kind=8), intent(in) :: c_hat(:)                                   !< Reacted mix
        !> Variables
            integer(kind=4) :: j                                                   !< Species index
            character(len=18) :: nm                                                !< Fixed-width (left-justified) row label so every row's value columns align under the species-name header
            logical, save :: descr_written = .false.  !< file variable-description block written once (first call)
        !> Process
            associate(sc => this%solid_chemistry)
            !> One-time file description: meaning of each labelled row within a per-cell block
            if (.not. descr_written) then
                write(iunit,'(A)') "[CELL-DBG] ==========================================================================="
                write(iunit,'(A)') "[CELL-DBG] REMIX per-cell reactive-mixing (WMA) trace  one block per (step, target water)"
                write(iunit,'(A)') "[CELL-DBG] concentrations [mol/L]; reaction amounts [mol/L per step]; rates [mol/L/time]"
                write(iunit,'(A)') "[CELL-DBG] variables per block:"
                write(iunit,'(A)') "[CELL-DBG]   aq var act species : variable-activity aqueous species names (column order)"
                write(iunit,'(A)') "[CELL-DBG]   components         : u_i as a linear combination of var-act species (step 1)"
                write(iunit,'(A)') "[CELL-DBG]   mineral eq reactions: mineral order for min re_mean / min amnt_eq (step 1)"
                write(iunit,'(A)') "[CELL-DBG]   c_mix   : transported (mixed) concentrations entering the reaction step"
                write(iunit,'(A)') "[CELL-DBG]   c_tilde : purely conservative mix (transport only, no reaction)"
                write(iunit,'(A)') "[CELL-DBG]   c_hat   : reacted mix = c_mix + Dt*r_tilde (= c_mix if lumped/no kinetics)"
                write(iunit,'(A)') "[CELL-DBG]   c_new   : concentrations after the equilibrium/kinetic solve (end of step)"
                write(iunit,'(A)') "[CELL-DBG]   u_mix / u_hat / u_new : component totals U*c for c_mix / c_hat / c_new"
                write(iunit,'(A)') "[CELL-DBG]   surf compl / surf act : surface complex names and Gaines-Thomas activities"
                write(iunit,'(A)') "[CELL-DBG]   CEC check : sum(z*c) adsorbed equivalents vs stored exchange capacity (CEC)"
                write(iunit,'(A)') "[CELL-DBG]   contrib species : column order of the contrib_* rows (aq species, then surf)"
                write(iunit,'(A)') "[CELL-DBG]   contrib_eq  : equilibrium reaction contribution = c_new - c_hat - kin term"
                write(iunit,'(A)') "[CELL-DBG]   contrib_kin : kinetic reaction contribution to the species balance"
                write(iunit,'(A)') "[CELL-DBG]   min/exch amnt_eq : equilibrium reaction amounts Re (minerals, then exchange)"
                write(iunit,'(A)') "[CELL-DBG]   aq r_eq / min r_eq : equilibrium reaction rates = amnt_eq / Dt"
                write(iunit,'(A)') "[CELL-DBG]   aq/min rk_mean : step-averaged kinetic reaction rates"
                write(iunit,'(A)') "[CELL-DBG]   aq/min rk_new, rk_old : end-of-step and previous-step kinetic rates"
                write(iunit,'(A)') "[CELL-DBG] ==========================================================================="
                descr_written = .true.
            end if
            write(iunit,'(A,I0,A,I0)') "[CELL-DBG] step=", k, " target water index=", twi
            nm='aq var act species'; write(iunit,'("[CELL-DBG]   ",A," =",*(A18))') nm, &
                (trim(this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(j)))%name), &
                 j=1,size(this%ind_var_act_species))
            if (k==1) call write_celldbg_legend(this, iunit)
            !> c_mix/c_hat are conc_nc-ordered: aqueous secondaries sit at offset num_prim_species
            !> (AFTER the free site / solid-gas primaries), not contiguously after the aqueous
            !> primaries. Print aqueous primaries then aqueous secondaries so the free site is not
            !> shown in place of an eliminated aqueous secondary (e.g. calcite's ca+2), keeping the
            !> row consistent with c_tilde/c_new. Non-exchange systems (num_prim==num_aq_prim) are
            !> unaffected: the second slice is then simply the contiguous secondaries.
            nm='c_mix'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, &
                [c_mix(1:sc%reactive_zone%speciation_alg%num_aq_prim_species), &
                 c_mix(sc%reactive_zone%speciation_alg%num_prim_species+1: &
                       sc%reactive_zone%speciation_alg%num_prim_species+ &
                       sc%reactive_zone%speciation_alg%num_aq_sec_var_act_species)]
            !> c_tilde is the aqueous-only conservative mix (size num_aq_var_act_species):
            !> print the whole array; u_tilde is not defined for it.
            nm='c_tilde'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, c_tilde
            nm='c_hat'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, &
                [c_hat(1:sc%reactive_zone%speciation_alg%num_aq_prim_species), &
                 c_hat(sc%reactive_zone%speciation_alg%num_prim_species+1: &
                       sc%reactive_zone%speciation_alg%num_prim_species+ &
                       sc%reactive_zone%speciation_alg%num_aq_sec_var_act_species)]
            nm='c_new'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%concentrations(this%ind_var_act_species)
            nm='u_mix'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%compute_u_mix(c_mix(1:n_v))
            nm='u_hat'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%compute_u_mix(c_hat(1:n_v))
            nm='u_new'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, &
                matmul(sc%reactive_zone%speciation_alg%comp_mat, this%get_conc_var_act())
            !> Surface (exchange) complexes are stored as the LAST num_surf_compl entries of the
            !> solid-chemistry activity/concentration arrays (after any minerals); print their
            !> names and activities just after u_new so the exchanger state is visible per cell.
            if (sc%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                block
                    integer(kind=4) :: nsc, off, nex, i, z
                    real(kind=8), allocatable :: act_sc(:)
                    real(kind=8) :: tot
                    nsc = sc%reactive_zone%cat_exch_zone%num_surf_compl
                    nex = sc%reactive_zone%cat_exch_zone%num_exch_cats
                    off = sc%num_solids - nsc
                    !> Gaines-Thomas activities are equivalent fractions. The equivalent of each
                    !> surface complex is z_i * c_i (free site x- has z=1), and the activities are
                    !> normalised by the total exchanger capacity ÃŽÂ£ z_j c_j so that the free-site
                    !> and adsorbed-cation activities always sum to 1.
                    allocate(act_sc(nsc))
                    act_sc(1) = sc%concentrations(off+1)   !< free site x- (valence 1)
                    do i=1,nex
                        z = sc%reactive_zone%chem_syst%aq_phase%aq_species( &
                            sc%reactive_zone%cat_exch_zone%exch_cat_indices(i))%valence
                        act_sc(1+i) = z * sc%concentrations(off+1+i)
                    end do
                    tot = sum(act_sc)
                    if (tot>0d0) act_sc = act_sc / tot
                    nm='surf compl'; write(iunit,'("[CELL-DBG]   ",A," =",*(A18))') nm, &
                        (trim(sc%reactive_zone%cat_exch_zone%surf_compl(j)%name), j=1,nsc)
                    nm='surf act'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, act_sc
                    !> Consistency check: the total adsorbed equivalents ÃŽÂ£ z_j c_j should equal the
                    !> stored cation exchange capacity CEC. Print both plus their relative mismatch.
                    write(iunit,'("[CELL-DBG]   CEC check: sum(z*c)=",ES18.10," CEC=",ES18.10," rel.diff=",ES18.10)') &
                        tot, sc%CEC, abs(tot-sc%CEC)/max(abs(sc%CEC),1d-30)
                end block
            end if
            !> contrib_eq/contrib_kin are stored in conc_nc order (aqueous primaries, solid/gas
            !> primaries incl the free site, aqueous secondaries, surface complexes). The header is
            !> printed as [aqueous var-act species (ind_var_act), then surface complexes], so for an
            !> exchange system the free site (at num_prim_species) and an aqueous secondary (e.g.
            !> ca+2, at num_prim_species+1) would otherwise be shown in each other's column. Remap the
            !> conc_nc-ordered contrib arrays onto that header order before printing.
            block
                integer(kind=4) :: naq, nsc, naqp, naqsec, nprim, j, s
                integer(kind=4), allocatable :: disp_order(:)
                naq = sc%reactive_zone%speciation_alg%num_aq_var_act_species
                nsc = sc%reactive_zone%cat_exch_zone%num_surf_compl
                nm='contrib species'; write(iunit,'("[CELL-DBG]   ",A," =",*(A18))') nm, &
                    (trim(this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(j)))%name), &
                     j=1,size(this%ind_var_act_species)), &
                    (trim(sc%reactive_zone%cat_exch_zone%surf_compl(j)%name), j=1,nsc)
                if (sc%reactive_zone%speciation_alg%flag_cat_exch) then
                    naqp   = sc%reactive_zone%speciation_alg%num_aq_prim_species
                    naqsec = sc%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                    nprim  = sc%reactive_zone%speciation_alg%num_prim_species
                    allocate(disp_order(naq+nsc))
                    do j=1,naqp;   disp_order(j)=j;            end do      !> aqueous primaries
                    do j=1,naqsec; disp_order(naqp+j)=nprim+j; end do      !> aqueous secondaries
                    disp_order(naq+1)=nprim                                !> free site (solid primary)
                    do s=2,nsc;    disp_order(naq+s)=nprim+naqsec+(s-1); end do !> adsorbed cations
                    nm='contrib_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%contrib_rxn_amnt_eq(disp_order)
                    if (allocated(this%contrib_rxn_amnt_kin)) then
                        nm='contrib_kin'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%contrib_rxn_amnt_kin(disp_order)
                    end if
                    deallocate(disp_order)
                else
                    nm='contrib_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%contrib_rxn_amnt_eq
                    if (allocated(this%contrib_rxn_amnt_kin)) then
                        nm='contrib_kin'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%contrib_rxn_amnt_kin
                    end if
                end if
            end block
            !> sc%rxn_amnt_eq and sc%rxn_rate_eq are combined arrays: minerals occupy entries
            !> 1..num_minerals and exchange reactions entries num_minerals+1..+num_exch_cats.
            !> Print them as separate rows so the mineral and exchange contributions are distinct.
            block
                integer(kind=4) :: nme, nex
                nme = sc%reactive_zone%num_minerals
                nex = sc%reactive_zone%cat_exch_zone%num_exch_cats
                if (nme>0) then
                    nm='min amnt_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_amnt_eq(1:nme)
                end if
                if (nex>0) then
                    nm='exch amnt_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_amnt_eq(nme+1:nme+nex)
                end if
                nm='aq r_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%rxn_rate_eq
                if (nme>0) then
                    nm='min r_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_rate_eq(1:nme)
                end if
                if (nex>0) then
                    nm='exch r_eq'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_rate_eq(nme+1:nme+nex)
                end if
            end block
            if (allocated(this%rxn_rate_kin_avg)) then
                nm='aq rk_mean'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%rxn_rate_kin_avg
            end if
            if (allocated(sc%rxn_rate_kin_avg)) then
                nm='min rk_mean'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_rate_kin_avg
            end if
            if (allocated(this%rxn_rate_kin)) then
                nm='aq rk_new'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%rxn_rate_kin
            end if
            if (allocated(this%rxn_rate_kin_old)) then
                nm='aq rk_old'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, this%rxn_rate_kin_old
            end if
            if (allocated(sc%rxn_rate_kin)) then
                nm='min rk_new'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_rate_kin
            end if
            if (allocated(sc%rxn_rate_kin_old)) then
                nm='min rk_old'; write(iunit,'("[CELL-DBG]   ",A," =",*(ES18.10))') nm, sc%rxn_rate_kin_old
            end if
            !> Flush so per-cell debug output survives even if a later time step aborts
            !> (e.g. convergence failure); otherwise buffered content is lost on error stop.
            flush(iunit)
            end associate
        end subroutine

!> Compute residual of component mass balance including constant activity species
!!> Mathematical Context: Component mass balance: U*c = u where U includes constant activity species
!!> The residual measures how well current concentrations satisfy component conservation
!!> Formula: residual = U*c - u (should be zero at equilibrium)
subroutine compute_res_spec_cst_act(this,conc_comp,conc,residual)
    implicit none                                                           !< No implicit variable declarations
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this                         !< Aqueous chemistry object (input only)
    real(kind=8), intent(in) :: conc_comp(:)                               !< Component concentrations (u vector)
    real(kind=8), intent(in) :: conc(:)                                    !< All species concentrations (c vector, includes constant activity)
    real(kind=8), intent(out) :: residual(:)                               !< Residual of component mass balance (pre-allocated)
    
    !> Process
    !> Compute component mass balance residual including constant activity species
    residual=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act,conc) - conc_comp
                                                                            !< residual = U_cst_act * c_all - u
end subroutine

!> Compute residual of component mass balance excluding constant activity species
!!> Mathematical Context: Component mass balance for variable activity species only
!!> This version uses only variable activity species in the mass balance calculation
!!> Formula: residual = U*c_nc - u where c_nc contains only variable activity species
subroutine compute_res_spec(this,conc_comp,c_nc,residual)
    implicit none                                                           !< No implicit variable declarations
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this                         !< Aqueous chemistry object (input only)
    real(kind=8), intent(in) :: conc_comp(:)                               !< Component concentrations (u vector)
    real(kind=8), intent(in) :: c_nc(:)                                    !< Variable activity species concentrations only
    real(kind=8), intent(out) :: residual(:)                               !< Residual of component mass balance (pre-allocated)
    
    !> Process
    !> Compute component mass balance residual for variable activity species only
    residual=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_nc) - conc_comp
                                                                            !< residual = U * c_nc - u
end subroutine

!> Updates concentration aqueous primary species in iterative method
subroutine update_conc_aq_prim_species(this,c1_aq,Delta_c1_aq,zero_flag)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: c1_aq(:) !> must be already allocated
    real(kind=8), intent(inout) :: Delta_c1_aq(:) !> must be already allocated
    logical, intent(out) :: zero_flag !> flag to indicate if concentration of aqueous primary species is zero

    integer(kind=4) :: i,n_p_aq
    real(kind=8), allocatable :: c1_aq_old(:)
    logical, save :: warned_zero_aq_prim = .false. !> emit the zero-concentration warning only once per run (avoids stdout flooding)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    zero_flag=.false. !> default value

    if (n_p_aq/=size(Delta_c1_aq) .and. this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species/=size(Delta_c1_aq)) then
        error stop "Dimension error in update_conc_aq_prim_species"
    end if
    c1_aq_old=c1_aq
    do i=1,n_p_aq
        if (c1_aq_old(i)+Delta_c1_aq(i)<=&
            this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            c1_aq_old(i)) then
            c1_aq(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
                c1_aq_old(i)
        else if (c1_aq_old(i)+Delta_c1_aq(i)>=c1_aq_old(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            c1_aq(i)=c1_aq_old(i)/&
                this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            c1_aq(i)=c1_aq_old(i)+Delta_c1_aq(i)
        end if
        if (c1_aq(i)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then
            !c1_aq(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol**2 !> chapuza
            if (.not. warned_zero_aq_prim) then
                print *, "Warning: concentration of aqueous primary species is zero (further warnings suppressed)"
                print *, "Aqueous primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_prim_species(i)))%name
                warned_zero_aq_prim=.true.
            end if
            zero_flag=.true. !> flag to indicate that concentration of aqueous primary species is zero
        end if
        Delta_c1_aq(i)=c1_aq(i)-c1_aq_old(i)
    end do
    call this%set_conc_aq_prim_species(c1_aq)
end subroutine

!> Updates concentration variable activity species in iterative method
subroutine update_conc_nc(this,conc_nc,Delta_c,zero_flag,sp_ind)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: conc_nc(:) !> concentration of variable activity species
    real(kind=8), intent(inout) :: Delta_c(:) !> difference in concentration of variable activity species
    logical, intent(out) :: zero_flag !> flag to indicate if concentration of variable activity species is zero
    integer(kind=4), intent(out) :: sp_ind !> index of species with zero concentration
    
    integer(kind=4) :: i,n_v,n_p_aq
    real(kind=8), allocatable :: c_old(:)
    logical, save :: warned_zero_prim = .false. !> emit the zero-concentration warning only once per run (avoids stdout flooding)
    
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !> number of variable activity species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> number of aqueus primary species
    zero_flag=.false. !> default value
    
    if (n_v/=size(Delta_c) .or. n_v/=size(conc_nc)) then
        error stop "Dimension error in update_conc_nc"
    end if
    c_old=conc_nc
    do i=1,n_v
        if (conc_nc(i)+Delta_c(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*conc_nc(i)) then
            conc_nc(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*conc_nc(i)
        else if (conc_nc(i)+Delta_c(i)>=conc_nc(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            conc_nc(i)=conc_nc(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            conc_nc(i)=conc_nc(i)+Delta_c(i)
        end if
        if (i<=n_p_aq .and. conc_nc(i)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then
            !conc_nc(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol**2 !> chapuza
            if (.not. warned_zero_prim) then
                print *, "Warning: concentration of primary species is zero (further warnings suppressed)"
                print *, "Primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(i)))%name
                warned_zero_prim=.true.
            end if
            zero_flag=.true. !> flag to indicate that concentration of variable activity species is zero
            sp_ind=i !> index of species with zero concentration
        end if
        Delta_c(i)=conc_nc(i)-c_old(i)
    end do
    call this%set_conc_var_act_species(conc_nc)
end subroutine

!> \brief Updates aqueous primary species concentrations with damping/control factor
!!
!! \details Applies a bounded update to aqueous primary species concentrations
!! using a control factor to prevent excessive changes in the Newton iteration.
!!
!! @param Delta_c1 concentration increment vector (modified on output to reflect actual change)
!!
subroutine update_conc_aq_prim_species_bis(this,Delta_c1)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq
    real(kind=8), allocatable :: c1_old(:)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    
    if (n_p_aq/=size(Delta_c1) .and. this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species/=size(Delta_c1)) then
        error stop "Dimension error in update_conc_aq_prim_species"
    end if
    c1_old=this%concentrations(1:n_p_aq)
    do i=1,n_p_aq
        if (this%concentrations(i)+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(i)) then
            this%concentrations(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(i)
        else if (this%concentrations(i)+Delta_c1(i)>=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(i)=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c1(i)
        end if
        Delta_c1(i)=this%concentrations(i)-c1_old(i)
    end do
end subroutine

!> Updates concentration aqueous species in Newton method
subroutine update_conc_aq_species(this,Delta_c_aq)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_aq(:)
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (this%solid_chemistry%reactive_zone%CV_params%control_factor>1d0 .or. &
    this%solid_chemistry%reactive_zone%CV_params%control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,this%aq_phase%num_species
        if (this%concentrations(i)+Delta_c_aq(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
        this%concentrations(i)) then
            this%concentrations(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_aq(i)>=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(i)=this%concentrations(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_aq(i)
        end if
        Delta_c_aq(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine update_conc_aq_species

!> Updates concentration aqueous and solid primary species in iterative method
subroutine update_conc_prim_species(this,c1,Delta_c1)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: c1(:)
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq,n_p
    real(kind=8), allocatable :: c1_old(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    
    if (n_p/=size(c1)) error stop "Dimension error in update_conc_prim_species"
        
    c1_old=c1 !> old primary concentrations
    do i=1,n_p
        if (c1_old(i)+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)) then
            c1(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)
        else if (c1_old(i)+Delta_c1(i)>=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            c1(i)=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            c1(i)=c1_old(i)+Delta_c1(i)
        end if
        Delta_c1(i)=c1(i)-c1_old(i)
    end do
    call this%set_conc_prim_species(c1)
end subroutine

!> This subroutine computes the concentration of variable activity species after mixing
!! we assume all mixing waters have the same aqueous species
subroutine compute_c_mix(this,conc_old_mix_waters,indices_aq_species_mix_waters,&
    mixing_ratios_conc,c_mix)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: this !> first mixing water (structural reference)
    real(kind=8), intent(in) :: conc_old_mix_waters(:,:) !> conc_old from each mixing water (num_species x num_mix_waters)
    integer(kind=4), intent(in) :: indices_aq_species_mix_waters(:,:) !> indices_aq_species from each mixing water (num_aq_species x num_mix_waters)
    real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as mix waters
    real(kind=8), intent(out), allocatable :: c_mix(:) !> concentration of variable activity species after mixing
!> Variables
    integer(kind=4) :: i, j, num_perm, num_mix_wat
    integer(kind=4) :: idx_aq, n_vas, n_aq_vas, n_aq_prim, n_aq_sec, n_prim
    logical :: found
!> Cache deep dereference chain
    associate(spec_alg => this%solid_chemistry%reactive_zone%speciation_alg, &
              ind_diss => this%aq_phase%ind_diss_solids, &
              ind_vas  => this%ind_var_act_species, &
              ind_aq   => this%indices_aq_species)
    n_vas     = spec_alg%num_var_act_species
    n_aq_vas  = spec_alg%num_aq_var_act_species
    n_aq_prim = spec_alg%num_aq_prim_species
    n_aq_sec  = spec_alg%num_aq_sec_var_act_species
    n_prim    = spec_alg%num_prim_species
    num_mix_wat = size(conc_old_mix_waters, 2)
!> Use block for stack-allocated automatic arrays (avoid heap allocation per call)
    block
    integer(kind=4) :: perm(n_vas, 2)
    real(kind=8) :: c_mix_aq(n_aq_vas)
    real(kind=8) :: c_mix_aux(n_aq_vas)
!> Initialise c_mix from old concentrations
    c_mix = this%get_conc_old_var_act_new_rz()
    num_perm = 0
    do i = 1, n_aq_vas
        idx_aq = ind_aq(ind_diss(i))
        !> First water contribution from "this"
        c_mix_aq(i) = mixing_ratios_conc(1) * this%conc_old(idx_aq)
        !> Build permutation table where ind_var_act_species differs from indices_aq_species
        if (ind_vas(i) /= idx_aq) then
            !> Check if i is already recorded in perm(:,2) Ã¢â‚¬â€ inline search
            found = .false.
            do j = 1, num_perm
                if (perm(j, 2) == i) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                num_perm = num_perm + 1
                perm(num_perm, 1) = i
                !> Inline get_ind_int: find idx_aq in ind_vas
                do j = 1, n_vas
                    if (ind_vas(j) == idx_aq) then
                        perm(num_perm, 2) = j
                        exit
                    end if
                end do
            end if
        end if
        !> Remaining waters contribution
        do j = 1, num_mix_wat
            c_mix_aq(i) = c_mix_aq(i) + mixing_ratios_conc(j+1) * &
                conc_old_mix_waters(indices_aq_species_mix_waters(ind_diss(i), j), j)
        end do
    end do
    !> Apply permutations via copy-swap
    c_mix_aux(1:n_aq_vas) = c_mix_aq(1:n_aq_vas)
    do i = 1, num_perm
        c_mix_aq(perm(i, 1)) = c_mix_aux(perm(i, 2))
        c_mix_aq(perm(i, 2)) = c_mix_aux(perm(i, 1))
    end do
    !> Copy results into c_mix using array slicing
    c_mix(1:n_aq_prim) = c_mix_aq(1:n_aq_prim)
    c_mix(n_prim+1:n_prim+n_aq_sec) = c_mix_aq(n_aq_prim+1:n_aq_prim+n_aq_sec)
    end block
    end associate
end subroutine

!> Consistent WMA global mixer: c_mix = conservative mix of all waters (aqueous var-act species) and
!> c_hat = c_mix in the full var-act layout plus the previous-step (1-theta_r) reaction term of the
!> mixing waters,  c_hat += omega*Delta_t*(1-theta_r)*sum_j lambda_j S_j,old^T r_old_j, with the uniform
!> positivity limiter omega in [0,1] of add_limited_increment
!> (react_lam is aligned with conc_idx: 0 for the waters that mix concentrations but not reaction rates).
subroutine compute_c_mix_global_cons(all_waters,zone_waters,conc_idx,conc_lam,react_lam,&
    wat_indices_zones,ext_wat_indices,Delta_t,theta_r,c_hat,c_mix,ideal)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:)
    type(zone_water_c), intent(inout) :: zone_waters(:) !> one temporary speciated water per zone
    integer(kind=4), intent(in) :: conc_idx(:)
    real(kind=8), intent(in) :: conc_lam(:), react_lam(:)
    type(int_array_c), intent(in) :: wat_indices_zones !> unused; kept for p_mixer interface compatibility
    integer(kind=4), intent(in) :: ext_wat_indices(:) !> unused; kept for p_mixer interface compatibility
    real(kind=8), intent(in) :: Delta_t, theta_r
    real(kind=8), intent(out), allocatable :: c_hat(:), c_mix(:)
    logical, intent(in) :: ideal !> unused (no speciation here); kept for p_mixer interface compatibility

    real(kind=8), allocatable :: c_mix_full(:) !> conservative mix in the full var-act layout (base of the limited increment)
    real(kind=8), allocatable :: r_tilde(:)

    call compute_c_mix_cons_global(all_waters,conc_idx,conc_lam,c_mix)
    associate(this => all_waters(conc_idx(1)))
    c_mix_full=expand_c_mix_cons(this,c_mix)
    c_hat=c_mix_full
    if (theta_r/=1d0) then
        allocate(r_tilde(size(c_hat)))
        r_tilde=0d0
        call add_prev_step_kin_react_mix_contrib(all_waters,conc_idx,size(conc_idx),react_lam,theta_r,r_tilde)
        !> Positivity-preserving update (same uniform-omega limiter as the other mixers)
        call add_limited_increment(c_mix_full,Delta_t*r_tilde,c_hat)
    end if
    end associate
end subroutine compute_c_mix_global_cons

!> Lumped WMA global mixer: c_mix = conservative mix of all waters (aqueous var-act species) and c_hat =
!> c_mix in the full var-act layout plus the SELF water's previous-step kinetic term weighted by the
!> lumped reaction mixing ratio,  c_hat += omega*lambda_L*Delta_t*(1-theta_r)*Sk_i,old^T rk_i,old,
!> with the uniform positivity limiter omega in [0,1] of add_limited_increment.
subroutine compute_c_mix_global_lump(all_waters,zone_waters,conc_idx,conc_lam,lambda_L,&
    wat_indices_zones,ext_wat_indices,Delta_t,theta_r,c_hat,c_mix,ideal)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:)
    type(zone_water_c), intent(inout) :: zone_waters(:) !> one temporary speciated water per zone
    integer(kind=4), intent(in) :: conc_idx(:)
    real(kind=8), intent(in) :: conc_lam(:)
    real(kind=8), intent(in) :: lambda_L !> lumped reaction mixing ratio of the target water
    type(int_array_c), intent(in) :: wat_indices_zones !> unused; kept for p_mixer interface compatibility
    integer(kind=4), intent(in) :: ext_wat_indices(:) !> unused; kept for p_mixer interface compatibility
    real(kind=8), intent(in) :: Delta_t, theta_r
    real(kind=8), intent(out), allocatable :: c_hat(:), c_mix(:)
    logical, intent(in) :: ideal !> unused (no speciation here); kept for p_mixer interface compatibility

    real(kind=8), allocatable :: c_mix_full(:) !> conservative mix in the full var-act layout (base of the limited increment)
    real(kind=8), allocatable :: Sk_nc_old(:,:), rk_old(:)
    integer(kind=4) :: p !> [NAN-DBG] row counter
    logical, save :: mix_dbg_done=.false. !> [MIX-DBG] one-shot flag of this call site

    if (within_zone_dbg_step==1 .and. .not. mix_dbg_done) then
        mix_dbg_done=.true.
        call print_mix_dbg_stencil("compute_c_mix_global_lump",conc_idx,conc_lam,lambda_L)
    end if
    call compute_c_mix_cons_global(all_waters,conc_idx,conc_lam,c_mix)
    associate(this => all_waters(conc_idx(1)))
    c_mix_full=expand_c_mix_cons(this,c_mix)
    c_hat=c_mix_full
    if (theta_r/=1d0 .and. this%indices_rk_old%num_cols>0) then
        Sk_nc_old=this%get_Sk_nc_old()
        rk_old=this%get_rk_old()
        !> Positivity-preserving update (same uniform-omega limiter as the other mixers)
        if (size(Sk_nc_old,2)==size(c_hat)) &
            call add_limited_increment(c_mix_full,Delta_t*(1d0-theta_r)*lambda_L*matmul(transpose(Sk_nc_old),rk_old),c_hat)
        !> [NAN-DBG] non-finite reacted mix: dump the kinetic previous-step term ingredients
        if (any(c_hat/=c_hat)) then
            print '(A,I0,A,I0,A,I0,A,I0,A,I0)', "[NAN-DBG] compute_c_mix_global_lump: target water ", conc_idx(1), &
                " indices_rk%num_cols=", this%indices_rk%num_cols, " indices_rk_old%num_cols=", this%indices_rk_old%num_cols, &
                " size(rk_old)=", size(rk_old), " rows(Sk_nc_old)=", size(Sk_nc_old,1)
            print '(A,*(ES12.4))', "[NAN-DBG]   c_mix_full =", c_mix_full
            print '(A,*(ES12.4))', "[NAN-DBG]   rk_old     =", rk_old
            do p=1,size(Sk_nc_old,1)
                print '(A,I0,A,*(ES12.4))', "[NAN-DBG]   Sk_nc_old row ", p, " =", Sk_nc_old(p,:)
            end do
            print '(A,*(ES12.4))', "[NAN-DBG]   c_hat      =", c_hat
            error stop "[NAN-DBG] non-finite reacted mix in compute_c_mix_global_lump"
        end if
    end if
    end associate
end subroutine compute_c_mix_global_lump

!> Expands the aqueous-only [aq prim | aq sec var act] conservative mix c_mix into the full
!> variable-activity layout of get_conc_old_var_act_new_rz; non-aqueous slots keep the water's previous-step values.
function expand_c_mix_cons(this,c_mix) result(c_hat)
    implicit none
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c_mix(:)
    real(kind=8), allocatable :: c_hat(:)

    integer(kind=4) :: n_aq_prim, n_prim, n_aq_sec

    n_aq_prim = this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_prim    = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_aq_sec  = this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
    c_hat = this%get_conc_old_var_act_new_rz()
    c_hat(1:n_aq_prim) = c_mix(1:n_aq_prim)
    c_hat(n_prim+1:n_prim+n_aq_sec) = c_mix(n_aq_prim+1:n_aq_prim+n_aq_sec)
end function expand_c_mix_cons

!> Conservative mixing of all waters: aqueous variable-activity species only, ordered
!> [aqueous primary | aqueous secondary var act]; also stored on the target water.
subroutine compute_c_mix_cons_global(all_waters,conc_idx,conc_lam,c_mix)
    implicit none
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:)
    integer(kind=4), intent(in) :: conc_idx(:) !> concentration-mixing waters (self first), global indices
    real(kind=8), intent(in) :: conc_lam(:) !> concentration mixing ratios aligned with conc_idx
    real(kind=8), intent(out), allocatable :: c_mix(:)

    integer(kind=4) :: species_pos, mix_pos, water_idx, aq_idx
    integer(kind=4) :: n_aq, n_aq_prim, n_aq_sec

    associate(this => all_waters(conc_idx(1)))
    n_aq      = this%aq_phase%num_species
    n_aq_prim = this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_aq_sec  = this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species

    block
    real(kind=8) :: c_tmp(n_aq)
    c_tmp=0d0
    do species_pos = 1, n_aq
        aq_idx = this%indices_aq_phase(species_pos)
        do mix_pos = 1, size(conc_idx)
            water_idx = conc_idx(mix_pos)
            c_tmp(species_pos) = c_tmp(species_pos) + conc_lam(mix_pos) * &
                all_waters(water_idx)%conc_old(all_waters(water_idx)%indices_aq_species(aq_idx))
        end do
    end do
    !> [NAN-DBG] non-finite conservative mix: report the stencil waters whose previous-step
    !> concentrations are non-finite (uninitialised / corrupted conc_old), then stop
    if (any(c_tmp /= c_tmp)) then
        print '(A,I0,A,*(ES12.4))', "[NAN-DBG] compute_c_mix_cons_global: target water ", conc_idx(1), &
            " non-finite mix, c_tmp (aq_phase order)=", c_tmp
        do mix_pos = 1, size(conc_idx)
            water_idx = conc_idx(mix_pos)
            if (any(all_waters(water_idx)%conc_old /= all_waters(water_idx)%conc_old) .or. conc_lam(mix_pos) /= conc_lam(mix_pos)) then
                print '(A,I0,A,I0,A,ES12.4,A,*(ES12.4))', "[NAN-DBG]   water ", water_idx, " (id ", all_waters(water_idx)%id, &
                    ") lambda_c=", conc_lam(mix_pos), " conc_old=", all_waters(water_idx)%conc_old
                print '(A,*(ES12.4))', "[NAN-DBG]     concentrations=", all_waters(water_idx)%concentrations
            end if
        end do
        error stop "[NAN-DBG] non-finite conservative mix"
    end if
    if (allocated(c_mix)) deallocate(c_mix)
    allocate(c_mix(n_aq_prim+n_aq_sec))
    do species_pos = 1, n_aq_prim
        c_mix(species_pos) = c_tmp(this%ind_prim_species(species_pos))
    end do
    do species_pos = 1, n_aq_sec
        c_mix(n_aq_prim + species_pos) = c_tmp(this%ind_sec_species(species_pos))
    end do
    end block
    this%conc_mix=c_mix
    end associate
end subroutine compute_c_mix_cons_global

!> This subroutine computes the concentration of aqueous variable activity species & aqueous kinetic reaction rates after water mixing
!! We assume all kinetic reactions are aqueous
!! we assume all mixing waters have the same kinetic reactions
!> COMMENTED OUT (dead code, no callers): compute_c_and_rk_tilde
!subroutine compute_c_and_rk_tilde(mix_waters,mixing_ratios_conc,mixing_ratios_R,c_mix,rk_tilde)
!    implicit none
!!> Argumnets
!    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
!    real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
!    real(kind=8), intent(in) :: mixing_ratios_R(:)  !> same order as "mix_waters"
!    real(kind=8), intent(out) :: c_mix(:) !> concentration of aqueous variable activity species after mixing
!    real(kind=8), intent(out) :: rk_tilde(:) !> aqueous kinetic reaction rates after mixing
!!> Variables
!    integer(kind=4) :: i !> index in "ind_var_act_species"
!    integer(kind=4) :: j !> index mixing waters
!    integer(kind=4) :: k !> index mineral kinetic reactions
!    integer(kind=4) :: num_perm !> counter in "perm"
!    integer(kind=4) :: num_mix_wat !> number of mixing waters
!    integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
!    real(kind=8), allocatable :: c_mix_aux(:) !> auxiliary variable to store the concentration of aqueous variable activity species after mixing
!    logical :: flag !> flag to check if an index is in "perm"
!!> Pre-process
!    allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
!    num_mix_wat=size(mix_waters)
!!> Process
!    c_mix=0d0 !> initialise concentration of aqueous variable activity species after mixing
!    rk_tilde=0d0 !> initialise aqueous kinetic reaction rates after mixing
!    num_perm=0 !< counter in "perm"
!    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
!        !> Autentica chapuza
!        if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
!            if (num_perm>0) then
!                call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
!                if (flag .eqv. .true.) then
!                    continue
!                else
!                    num_perm=num_perm+1
!                    perm(num_perm,1)=i
!                    perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
!                        mix_waters(1)%ind_var_act_species)
!                end if
!            else
!                num_perm=num_perm+1
!                perm(num_perm,1)=i
!                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
!                    mix_waters(1)%ind_var_act_species)
!            end if
!        end if
!        do j=1,num_mix_wat
!            c_mix(i)=c_mix(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
!                mix_waters(j)%aq_phase%ind_diss_solids(i)))
!
!        end do
!    end do
!    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
!        do j=1,num_mix_wat
!            rk_tilde(i)=rk_tilde(i)+mixing_ratios_R(j)*mix_waters(j)%rxn_rate_kin(i)
!        end do
!    end do
!    if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
!            j=1
!            k=1
!            i=1
!            do
!                if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
!                    mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
!                    rk_tilde(i)=rk_tilde(i)+mixing_ratios_R(j)*mix_waters(j)%solid_chemistry%rxn_rate_kin(i)
!                    if (j<num_mix_wat) then
!                        j=j+1
!                        i=1
!                        k=1
!                    else
!                        exit
!                    end if
!                else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
!                    k=k+1
!                else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
!                    i=i+1
!                    k=1
!                else if (j<num_mix_wat) then
!                    j=j+1
!                    i=1
!                    k=1
!                else
!                    exit
!                end if
!            end do
!    end if
!    c_mix_aux=c_mix !> we store the concentration of aqueous variable activity species after mixing (chapuza)
!    do i=1,num_perm
!        !> we permute the indices in "perm"
!        c_mix(perm(i,1))=c_mix_aux(perm(i,2))
!        c_mix(perm(i,2))=c_mix_aux(perm(i,1))
!    end do
!end subroutine 


!> This subroutine compares the associated reactions of two aqueous chemistry objects.
!> Two waters belong to the same reaction zone when they share the same equilibrium reactions
!> (identified by the non-flowing species of their reactive zone) and the same kinetic reactions
!> (identified by the kinetic minerals of their mineral zone). Aqueous kinetic reactions are shared
!> through the chemical system, so they do not distinguish reaction zones.
subroutine compare_old_reaction_zones(aq_chem_1,aq_chem_2,flag)
    use reactive_zone_m, only: compare_react_zones
    implicit none
    class(aqueous_chemistry_c), intent(in) :: aq_chem_1 !> first water
    class(aqueous_chemistry_c), intent(in) :: aq_chem_2 !> second water
    logical, intent(out) :: flag !> TRUE if same reaction zone (equilibrium + kinetic reactions), FALSE otherwise

    integer(kind=4) :: n_kin_1 !> number of kinetic minerals of the first water
    integer(kind=4) :: n_kin_2 !> number of kinetic minerals of the second water

    !> Compare equilibrium reactions (non-flowing species of the reactive zone)
    call compare_react_zones(aq_chem_1%solid_chemistry_old%reactive_zone, &
                             aq_chem_2%solid_chemistry_old%reactive_zone,flag)
    if (.not. flag) return !> different equilibrium reactions => different reaction zone

    !> Number of kinetic minerals of each water (0 if no mineral zone is associated)
    n_kin_1=0
    if (associated(aq_chem_1%solid_chemistry_old%mineral_zone)) &
        n_kin_1=aq_chem_1%solid_chemistry_old%mineral_zone%num_minerals_kin
    n_kin_2=0
    if (associated(aq_chem_2%solid_chemistry_old%mineral_zone)) &
        n_kin_2=aq_chem_2%solid_chemistry_old%mineral_zone%num_minerals_kin

    !> Compare kinetic reactions (kinetic minerals of the mineral zone)
    if (n_kin_1/=n_kin_2) then
        flag=.false.
    else if (n_kin_1>0) then
        flag=all(aq_chem_1%solid_chemistry_old%mineral_zone%ind_min_chem_syst(1:n_kin_1)== &
                 aq_chem_2%solid_chemistry_old%mineral_zone%ind_min_chem_syst(1:n_kin_1))
    end if
end subroutine

!> \brief Adds the previous-time-step (1-theta) reactive-mixing contribution to r_tilde, summed over
!>        the reaction-mixing waters, the current water included (NO reactive-zone grouping).
!> \details   r_tilde += sum_m lambda_m [ Se_m,old^T re_old_m + (1-theta) Sk_m,old^T rk_old_m ],
!>          each water m (self and neighbours) with its OWN previous-step FULL stoichiometry S_old
!>          (EQUILIBRIUM rows, then KINETIC rows: get_S_nc_old) and old reaction rates r_old (same row
!>          order: get_r_old), weighted by its reaction mixing ratio. Only the KINETIC rates carry the
!>          time weighting (1-theta); the equilibrium rates enter with full weight (the equilibrium
!>          problem is not time-weighted). The equilibrium block is included because, without zone
!>          grouping, the mixing waters need not share the target's component matrix (U*S_eq^T = 0
!>          only holds within one reactive zone). The equilibrium columns are those of the TARGET's
!>          (self, mix_idx(1)) current reactive zone, the basis of r_tilde.
    subroutine add_prev_step_react_mix_contrib(all_waters,mix_idx,num_mix,mixing_ratios_reacts,&
        theta,r_tilde)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: mix_idx(:) !> global indices of the reaction-mixing waters, self first
    integer(kind=4), intent(in) :: num_mix !> number of reaction-mixing waters (self included)
    real(kind=8), intent(in) :: mixing_ratios_reacts(:) !> reaction mixing ratios, position-aligned with mix_idx
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: r_tilde(:) !> reaction rate contributions (updated in place)

    integer(kind=4) :: p !> mixing water position
    integer(kind=4) :: m !> global index of a mixing water
    integer(kind=4) :: n_eq_m !> number of equilibrium reactions of water m's previous-step zone (leading rows of S_m_old / r_old_m)
    real(kind=8) :: lambda_m !> reaction mixing ratio of mixing water m
    integer(kind=4), allocatable :: ind_var_act_species_rz(:) !> variable activity species indices of the target's current reactive zone (equilibrium columns of get_S_nc_old)
    real(kind=8), allocatable :: S_m_old(:,:) !> previous-step FULL stoichiometric matrix of a water (equilibrium rows, then kinetic)
    real(kind=8), allocatable :: r_old_m(:) !> previous-step FULL reaction rates of a water (equilibrium, then kinetic; get_r_old order)

    ind_var_act_species_rz=all_waters(mix_idx(1))%solid_chemistry%reactive_zone%ind_var_act_species
    do p=1,num_mix
        m=mix_idx(p)
        lambda_m=mixing_ratios_reacts(p)
        if (lambda_m==0d0) cycle
        S_m_old=all_waters(m)%get_S_nc_old(ind_var_act_species_rz)
        r_old_m=all_waters(m)%get_r_old()
        !> a water whose previous-step zone has a different variable-activity basis cannot be mapped
        !> onto r_tilde (same convention as compute_c_mix_by_zones_cons)
        if (size(S_m_old,2)/=size(r_tilde)) cycle
        !> time weighting (1-theta) on the KINETIC rows only; the equilibrium rows keep full weight
        !> (the kinetic tail of get_r_old has indices_rk_old%num_cols rows, see get_r_hist)
        n_eq_m=size(r_old_m)-all_waters(m)%indices_rk_old%num_cols
        r_old_m(n_eq_m+1:)=(1d0-theta)*r_old_m(n_eq_m+1:)
        r_tilde=r_tilde+lambda_m*matmul(transpose(S_m_old),r_old_m)
    end do
end subroutine

!> \brief Same as add_prev_step_react_mix_contrib but for mixing waters that all share the same
!>        reactive zone: the EQUILIBRIUM reactions are omitted (KINETIC reactions only).
!> \details When every mixing water belongs to the same reactive zone, the equilibrium part of the
!>          previous-step reaction-mixing term is annihilated by the shared component matrix
!>          (U*S_eq^T = 0), so it can be dropped; only each water's KINETIC reactions are summed:
!>            r_tilde += (1-theta) * sum_m lambda_m Sk_m,old^T r_kin,old_m.
!>          get_Sk_nc_old and get_rk_old both return the KINETIC block directly (no
!>          variable-activity-species index or equilibrium-block slicing needed).
    subroutine add_prev_step_kin_react_mix_contrib(all_waters,mix_idx,num_mix,&
        mixing_ratios_reacts,theta,r_tilde)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: mix_idx(:) !> global indices of the mixing waters
    integer(kind=4), intent(in) :: num_mix !> number of mixing waters
    real(kind=8), intent(in) :: mixing_ratios_reacts(:) !> reaction mixing ratios, position-aligned with mix_idx
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: r_tilde(:) !> reaction rate contributions (updated in place)

    integer(kind=4) :: p !> mixing water position
    integer(kind=4) :: m !> global index of a mixing water
    real(kind=8) :: lambda_m !> reaction mixing ratio of mixing water m
    real(kind=8), allocatable :: Sk_m_old(:,:) !> mixing water's previous-step KINETIC stoichiometric matrix
    real(kind=8), allocatable :: rk_old_m(:) !> mixing water's previous-step KINETIC reaction rates

    !> Previous-step (1-theta) contribution summed directly over the mixing waters, KINETIC
    !> reactions only (the equilibrium block cancels through the shared component matrix):
    !>   r_tilde += (1-theta) * sum_m lambda_m Sk_m,old^T rk_old_m
    do p=1,num_mix
        m=mix_idx(p)
        lambda_m=mixing_ratios_reacts(p)
        if (lambda_m==0d0) cycle
        rk_old_m=all_waters(m)%get_rk_old()
        if (size(rk_old_m)==0) cycle !> no kinetic reactions in this neighbour
        Sk_m_old=all_waters(m)%get_Sk_nc_old()
        r_tilde=r_tilde+(1d0-theta)*lambda_m*matmul(transpose(Sk_m_old),rk_old_m)
    end do
end subroutine

!> \brief Page-4 within-reactive-zone mixing (algorithm step 2), shared worker of
!>        mix_within_zone_old_rz / mix_within_zone_new_rz: produces a TEMPORARY water holding the
!>        speciated sub-mixture (the "mixed water by reactive zone") formed by the mixing waters of the
!>        target cell that are all grouped in the SAME reactive zone z.
!> \details Zone z is defined by rep_rz, a reactive zone of the group's representative (first member):
!>          its OLD one (solid_chemistry_old%reactive_zone: the zone the members' previous-step
!>          concentrations were equilibrated with; for a displaced Lagrangian parcel, the zone of the
!>          cell it comes from) or its CURRENT one (solid_chemistry%reactive_zone: the zone of the cell
!>          it sits in). The sub-mixture is speciated on a scratch solid chemistry copied from the
!>          representative's CURRENT solid chemistry (mineral zone, exchanger state) and re-pointed to
!>          rep_rz.
!>            c_mix_z = sum_m lambda_c,z,m c_m^k + Delta_t*(1-theta) sum_m lambda_r,z,m Sk_m^T rk_m^k (limited),
!>          where the previous-step reaction term is KINETIC only (add_prev_step_kin_react_mix_contrib):
!>          the equilibrium block cancels under the zone's component matrix (U_z*S_eq^T = 0) and is
!>          never added. c_mix_z is then speciated with zone z (equilibrate_water(_ideal) on a temporary
!>          copy of the representative moved to zone z, leaving the real water states untouched) so its
!>          aqueous species can later be re-mixed across zones by the target cell - but ONLY if rep_rz
!>          differs from the reactive zone the members are equilibrated with (the representative's own
!>          solid_chemistry%reactive_zone): a group already in rep_rz keeps its conservative mix, as
!>          the target's own equilibration conserves the same components. The within-zone
!>          ratios must already be renormalised to sum to 1 over the group
!>          (lambda_*,z,m = lambda_*,i,m / sum_z lambda_*,i,n).
subroutine mix_within_zone_core(this,all_waters,zone_mix_idx,ratios_conc_norm,ratios_reacts_norm,&
    Delta_t,theta_r,ideal,rep_rz,rep_mz)
    implicit none
    class(aqueous_chemistry_c), intent(inout) :: this !> water of zone z holding the speciated sub-mixture
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: zone_mix_idx(:) !> global indices of the group's mixing waters (all grouped in zone z)
    real(kind=8), intent(in) :: ratios_conc_norm(:) !> within-zone concentration mixing ratios (sum=1), aligned with zone_mix_idx
    real(kind=8), intent(in) :: ratios_reacts_norm(:) !> within-zone reaction mixing ratios (sum=1), aligned with zone_mix_idx
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta_r !> time weighting factor
    logical, intent(in) :: ideal !> .true. ideal speciation (mixing_iter_comp_ideal); .false. non-ideal (mixing_iter_comp)
    class(reactive_zone_c), intent(in), target :: rep_rz !> reactive zone of water where final speciation occurs; target: the scratch solid chemistry points to it
    class(mineral_zone_c), intent(in), target :: rep_mz !> mineral zone of the water where final speciation occurs
    !type(solid_chemistry_c), intent(inout) :: zone_solid !> copy of the sub-mixture's scratch solid chemistry (its eq rates), taken before it is discarded

    integer(kind=4) :: p !> group member position
    integer(kind=4) :: n_v_z !> number of variable activity species in zone z
    integer(kind=4) :: n_p_z !> number of primary species in zone z
    integer(kind=4) :: n_aq_p_z !> number of aqueous primary species in zone z
    integer(kind=4) :: n_aq_s_z !> number of aqueous secondary variable activity species in zone z
    real(kind=8), allocatable :: c_mix_z(:) !> conservative within-zone mix of the previous-step variable activity concentrations
    real(kind=8), allocatable :: r_tilde_z(:) !> within-zone kinetic previous-step term
    real(kind=8), allocatable :: c_hat_z(:) !> pre-speciation concentrations: c_mix_z + omega*Delta_t*r_tilde_z, omega in [0,1] from add_limited_increment (speciation input)
    real(kind=8), allocatable :: c_eq_z(:) !> speciated (equilibrated) variable activity concentrations (Newton output; state stored in this)
    real(kind=8), allocatable :: conc_comp_z(:) !> component concentrations (work array for the speciation)
    real(kind=8), allocatable :: c1_old_init(:) !> primary-species initial guess for the Newton speciation
    logical :: speciate !> .true. if the sub-mixture is taken to a reactive zone other than its members' own (rep_rz /= the representative's zone): then it is speciated in rep_rz; otherwise it stays the conservative mix
    !class(solid_chemistry_c), pointer :: zone_solid_chem => null() !> temporary scratch solid chemistry for the zone sub-mixture (discarded on exit)

    !> The group's first water represents zone z. The members may sit in cells of OTHER zones (displaced
    !> parcels), so nothing zone-dependent is taken from their own solid chemistry or index arrays.
    associate(rep_water => all_waters(zone_mix_idx(1)), zone_solid => this%solid_chemistry)
    n_v_z=rep_rz%speciation_alg%num_var_act_species
    n_p_z=rep_rz%speciation_alg%num_prim_species
    n_aq_p_z=rep_rz%speciation_alg%num_aq_prim_species
    n_aq_s_z=rep_rz%speciation_alg%num_aq_sec_var_act_species
    !> Build the temporary zone water from the representative (speciated in place below) and move it
    !> to zone z: it gets its OWN scratch solid chemistry, a copy of the representative's CURRENT
    !> solid chemistry (mineral zone and exchanger state) re-pointed to rep_rz (zone z's reactive
    !> zone), so the speciation (speciation algebra, stoichiometry,
    !> logK) and compute_Re (which stores the sub-mixture's equilibrium rates, amounts and
    !> contributions on this and this scratch solid chemistry) all see zone z and never clobber a
    !> real solid chemistry. The index arrays are rebuilt in zone z's column order from the parcel's
    !> own aqueous-species indices, and the per-water logK / logK_star for zone z's reactions at the
    !> water's temperature.
    call this%copy_aqueous_attributes(rep_water)
    !allocate(zone_solid_chem)
    call zone_solid%set_reactive_zone(rep_rz)
    call zone_solid%set_mineral_zone(rep_mz)
    call zone_solid%initialise_solid_chemistry()
    !> The copied solid state (solid concentrations, eq rate/amount arrays: copy_local_chemistry
    !> carries them over, so the allocator above keeps their sizes) must be laid out for rep_rz; a
    !> representative sitting in another zone would otherwise make compute_Re / the speciation write
    !> past these arrays (silent heap corruption). Error stop instead.
    if (size(zone_solid%rxn_amnt_eq)/=rep_rz%num_minerals+rep_rz%cat_exch_zone%num_exch_cats .or. &
        zone_solid%num_solids/=rep_rz%num_solids+zone_solid%mineral_zone%num_minerals_kin) then
        print '(A,I0,A,I0,A,I0,A)', "[WITHIN-ZONE] ERROR: zone rep=", zone_mix_idx(1), &
            " its solid chemistry (zone id ", rep_water%solid_chemistry%reactive_zone%id, &
            ") is not consistent with rep_rz (zone id ", rep_rz%id, ")"
        error stop "[WITHIN-ZONE] representative's solid chemistry not consistent with the zone's reactive zone"
    end if
    !call this%set_solid_chemistry(zone_solid)
    !call this%set_solid_chemistry_old()
    call this%set_attributes_dep_rz()
    call this%set_attributes_dep_mz()
    !> Conservative within-zone mixing of the previous-step AQUEOUS variable activity concentrations,
    !> in zone z's (rep_rz's) column order: every member is read through conc_old_var_act_rep_rz,
    !> whatever the layout of the zone it sits in (a member's own ind_* arrays follow ITS current zone,
    !> whose primary/secondary split may differ from rep_rz's). The non-aqueous slots of c_mix_z
    !> (exchange free site, adsorbed cations, gases) stay 0.
    allocate(c_mix_z(n_v_z))
    c_mix_z=0d0
    do p=1,size(zone_mix_idx)
        c_mix_z=c_mix_z+ratios_conc_norm(p)*conc_old_var_act_rep_rz(all_waters(zone_mix_idx(p)))
    end do
    if (.not. allocated(this%conc_mix)) call this%allocate_conc_mix()
    this%conc_mix(1:n_aq_p_z)=c_mix_z(1:n_aq_p_z)
    this%conc_mix(n_aq_p_z+1:n_aq_p_z+n_aq_s_z)=c_mix_z(n_p_z+1:n_p_z+n_aq_s_z)
    !> Previous-step reaction term: kinetic only (the equilibrium block cancels under U_z)
    allocate(r_tilde_z(n_v_z))
    r_tilde_z=0d0
    call add_prev_step_kin_react_mix_contrib(all_waters,zone_mix_idx,size(zone_mix_idx),&
        ratios_reacts_norm,theta_r,r_tilde_z)
    allocate(c_hat_z(n_v_z))
    !> Positivity-preserving update (same uniform-omega limiter as the mixers)
    call add_limited_increment(c_mix_z,Delta_t*r_tilde_z,c_hat_z)
    !> the sub-mixture's own aqueous mix entering the speciation, c_hat_z ([aq prim | aq sec]: the
    !> conservative mix plus the limited previous-step kinetic term), for the output file and for
    !> the cross-zone c_mix of stage 3: two slice copies into the preallocated array (no
    !> array-constructor temporary)
    if (.not. allocated(this%conc_hat)) call this%allocate_conc_hat()
    this%conc_hat(1:n_aq_p_z)=c_hat_z(1:n_aq_p_z)
    this%conc_hat(n_aq_p_z+1:n_aq_p_z+n_aq_s_z)=c_hat_z(n_p_z+1:n_p_z+n_aq_s_z)
    allocate(c_eq_z(n_v_z))
    allocate(conc_comp_z(n_p_z))
    !> Newton initial guess: the primary block of the sub-mixture itself. c_hat_z is the conservative
    !> mix c_mix_z of equilibrated waters plus the (limited) previous-step kinetic term Delta_t*r_tilde_z,
    !> so it is NOT conservative; still, its primaries are close to the speciated ones (the kinetic
    !> increment is small over one step, and the within-zone mix only has to be re-equilibrated with
    !> the zone's equilibrium reactions) and a much better guess than the representative member's own
    !> state. equilibrate_water(_ideal) reads the guess from get_c1(), so it is written into
    !> this's primaries first; the representative's current primaries (read from this in
    !> zone z's primary order) only serve as the anchor of the retry extrapolation (mu > 0).
    c1_old_init=this%get_c1()
    call this%set_conc_prim_species(c_hat_z(1:n_p_z))
    !> The sub-mixture is speciated only when it is taken to a reactive zone OTHER than the one its
    !> members are equilibrated with (rep_rz /= the representative's own reactive zone): a group that
    !> already belongs to rep_rz (the target's own zone) only needs its conservative mix, since the
    !> target's own equilibration conserves exactly the components the within-zone speciation would
    !> (U_z*c_eq_z = U_z*c_hat_z) - the same redundancy that skips stage 2 altogether in the direct
    !> mix. Its zone water then carries c_hat_z itself (no equilibrium reaction: compute_Re below finds
    !> rho = c - c_hat = 0) and the equilibrium check does not apply to it.
    speciate=(rep_water%solid_chemistry%reactive_zone%id/=rep_rz%id)
    if (within_zone_dbg_active) then
        write(within_zone_dbg_unit,'(A,I0,A,I0,A,I0,A,L1)') "[ZONE-DBG]   zone rep water=", zone_mix_idx(1), &
            " n_members=", size(zone_mix_idx), " n_var_act=", n_v_z, " speciated=", speciate
        write(within_zone_dbg_unit,'("[ZONE-DBG]     zone var act =",*(A13))') &
            (trim(rep_rz%chem_syst%species(rep_rz%ind_var_act_species(p))%name), p=1,n_v_z)
        write(within_zone_dbg_unit,'("[ZONE-DBG]     ratios_conc  =",*(ES13.5))') ratios_conc_norm
        write(within_zone_dbg_unit,'("[ZONE-DBG]     r_tilde_z    =",*(ES13.5))') r_tilde_z
        write(within_zone_dbg_unit,'("[ZONE-DBG]     c_mix_z      =",*(ES13.5))') c_mix_z
        write(within_zone_dbg_unit,'("[ZONE-DBG]     c_hat_z      =",*(ES13.5))') c_hat_z
    end if
    if (.not. speciate) then
        c_eq_z=c_hat_z
        conc_comp_z=matmul(rep_rz%speciation_alg%comp_mat,c_hat_z)
        !> the zone water's aqueous state is the conservative sub-mixture itself (stage 3 reads it
        !> through concentrations); the non-aqueous slots of c_hat_z are 0 and are left alone
        this%concentrations(this%ind_var_act_species(1:n_aq_p_z))=c_hat_z(1:n_aq_p_z)
        this%concentrations(this%ind_var_act_species(n_aq_p_z+1:n_aq_p_z+n_aq_s_z))=c_hat_z(n_p_z+1:n_p_z+n_aq_s_z)
    else if (ideal) then
        call equilibrate_water_ideal(this,c1_old_init,c_hat_z,0d0,Delta_t,theta_r,c_eq_z,conc_comp_z)
    else
        call equilibrate_water(this,c1_old_init,c_hat_z,0d0,Delta_t,theta_r,c_eq_z,conc_comp_z)
    end if
    if (within_zone_dbg_active) write(within_zone_dbg_unit,'("[ZONE-DBG]     c_eq_z       =",*(ES13.5))') c_eq_z
    !> Component mass-balance check: the speciation must conserve the zone's components,
    !> U_z*c_eq_z == conc_comp_z (= U_z*c_hat_z). Error stop otherwise.
    block
        real(kind=8) :: mb_resid !> max |U_z*c_eq_z - conc_comp_z| over the zone's components
        mb_resid=maxval(abs(matmul(rep_rz%speciation_alg%comp_mat,c_eq_z)-conc_comp_z))
        if (mb_resid>max(rep_rz%CV_params%abs_tol,1d-12)) then
            print '(A,I0,A,ES12.4)', "[WITHIN-ZONE] ERROR: zone rep=", zone_mix_idx(1), &
                " component mass balance violated, max|U*c_eq-u|=", mb_resid
            error stop "[WITHIN-ZONE] zone sub-mixture component mass balance violated"
        end if
    end block
    !> Equilibrium reaction rates, amounts and contributions of the sub-mixture, stored on this
    !> and its scratch solid chemistry. lambda_r=1: full c_hat_z -> c_eq_z transfer; theta_r=0: no
    !> new-step kinetic term (its previous-step part is already folded into c_hat_z).
    call this%compute_Re(c_hat_z, Delta_t, 1d0, 0d0)
    !> Check that this is equilibrated with its (associated) reactive zone: for every
    !> equilibrium reaction the ion activity product prod_i a_i^S_ri must equal K (the Newton
    !> speciation in mixing_iter_comp_ideal/mixing_iter_comp enforces it). Error stop otherwise.
    !> Only meaningful when the sub-mixture was speciated here (a group of rep_rz's own zone is
    !> left as its conservative mix, which is not at equilibrium).
    call this%compute_activities_aq()
    call this%compute_pH()
    if (speciate .and. this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
        block
            real(kind=8), allocatable :: K_eq(:) !> equilibrium constants of the zone's reactions
            real(kind=8), allocatable :: a_col(:) !> activities in stoich_mat COLUMN (speciation) order
            real(kind=8) :: IAP !> ion activity product prod_i a_i^S_ri of one reaction
            real(kind=8) :: eq_resid !> max |K/IAP - 1| over the zone's equilibrium reactions
            integer(kind=4) :: r, s, n_p_aq, n_s_va, n_prim, n_sp
            n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            n_s_va=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            n_prim=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            n_sp=this%solid_chemistry%reactive_zone%speciation_alg%num_species
            !> activities are in aq_phase order, stoich_mat columns in speciation order: remap the
            !> aqueous activities into column order via the primary/secondary index arrays; every
            !> other column (minerals a=1, ideal water a=1) stays 1 so it drops out of the product.
            allocate(a_col(n_sp)); a_col=1d0
            do s=1,n_p_aq
                a_col(s)=this%activities(this%ind_prim_species(s))
            end do
            do s=1,n_s_va
                a_col(n_prim+s)=this%activities(this%ind_sec_species(s))
            end do
            !> this water's own K(T): the speciation used this%logK, not the 25 C eq_cst_ref
            K_eq=1d1**this%logK(1:this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
            eq_resid=0d0
            do r=1,size(K_eq)
                IAP=1d0
                do s=1,n_sp
                    IAP=IAP*a_col(s)**this%solid_chemistry%reactive_zone%stoich_mat(r,s)
                end do
                eq_resid=max(eq_resid,abs(K_eq(r)/IAP-1d0))
            end do
            !> Floor the tolerance above IAP-product floating-point roundoff (~1e-15): a converged
            !> water sits at machine epsilon, real non-equilibrium is O(1); 1e-12 (verified: cc_anh
            !> full run passes at this floor) cleanly separates them.
            if (eq_resid>max(this%solid_chemistry%reactive_zone%CV_params%abs_tol,1d-12)) then
                print '(A,I0,A,ES12.4)', "[WITHIN-ZONE] ERROR: zone rep=", zone_mix_idx(1), &
                    " NOT equilibrated with its reactive zone, max|K/IAP-1|=", eq_resid
                print '(A,4(1X,I0))', "[WITHIN-ZONE]   n_eq_react,n_aq_sec_va,n_prim,n_minerals=", &
                    size(K_eq), n_s_va, n_prim, this%solid_chemistry%reactive_zone%num_minerals
                do r=1,size(K_eq)
                    IAP=1d0
                    do s=1,n_sp
                        IAP=IAP*a_col(s)**this%solid_chemistry%reactive_zone%stoich_mat(r,s)
                    end do
                    print '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)', "[WITHIN-ZONE]   react r=", r, &
                        " K=", K_eq(r), " IAP=", IAP, " K/IAP=", K_eq(r)/IAP
                end do
                error stop "[WITHIN-ZONE] zone sub-mixture not equilibrated with its reactive zone"
            end if
        end block
    end if
    if (within_zone_dbg_active) then
        write(within_zone_dbg_unit,'("[ZONE-DBG]     rho_eq (aq)  =",*(ES13.5))') this%contrib_rxn_amnt_eq
        write(within_zone_dbg_unit,'("[ZONE-DBG]     Re (aq)      =",*(ES13.5))') this%rxn_amnt_eq
        write(within_zone_dbg_unit,'("[ZONE-DBG]     re_mean (aq) =",*(ES13.5))') this%rxn_rate_eq
        write(within_zone_dbg_unit,'("[ZONE-DBG]     Re (solid)   =",*(ES13.5))') zone_solid%rxn_amnt_eq
        write(within_zone_dbg_unit,'("[ZONE-DBG]     re_mean(sol) =",*(ES13.5))') zone_solid%rxn_rate_eq
    end if
    !zone_solid_out=zone_solid !> copy the sub-mixture's scratch solid chemistry to the output variable
    !deallocate(zone_solid) !> temporary: discard the zone's scratch solid chemistry
    !nullify(this%solid_chemistry) !> avoid a dangling alias to the freed scratch object
    end associate

    contains

    !> Previous-step concentrations of the AQUEOUS variable activity species of group member wm, in
    !> ZONE z's (rep_rz's) column order: the layout c_mix_z and the speciation of the sub-mixture use.
    !> Only the aqueous species are mixed within the zone: the result has rep_rz's full variable
    !> activity layout, with the aqueous primary block (1:n_aq_p_z) and the aqueous secondary block
    !> (n_p_z+1:n_p_z+n_aq_s_z) filled and every non-aqueous slot (exchange free site, adsorbed
    !> cations, gases) left at 0. A member's own ind_* arrays follow ITS current zone, whose
    !> primary/secondary split may differ from rep_rz's, so they are not used: each column is mapped
    !> through the aqueous-phase species index instead. this%index arrays (rebuilt for rep_rz
    !> by set_attributes_dep_rz) give the aqueous species of each column, and wm's own
    !> indices_aq_species locates it in wm%conc_old.
    function conc_old_var_act_rep_rz(wm) result(c)
        type(aqueous_chemistry_c), intent(in) :: wm !> group member
        real(kind=8) :: c(n_v_z) !> its previous-step aqueous variable activity concentrations in zone z's layout (non-aqueous slots 0)
        integer(kind=4) :: j !> zone column counter
        integer(kind=4) :: aq_idx !> aqueous-phase index of the species in a zone column
        c=0d0
        do j=1,n_aq_p_z
            aq_idx=this%indices_aq_phase(this%ind_var_act_species(j))
            c(j)=wm%conc_old(wm%indices_aq_species(aq_idx))
        end do
        do j=1,n_aq_s_z
            aq_idx=this%indices_aq_phase(this%ind_var_act_species(n_aq_p_z+j))
            c(n_p_z+j)=wm%conc_old(wm%indices_aq_species(aq_idx))
        end do
    end function conc_old_var_act_rep_rz
end subroutine mix_within_zone_core

!> \brief Within-zone mixing in the OLD reactive zone of the group (solid_chemistry_old%reactive_zone
!>        of its representative): the zone the members' previous-step concentrations were equilibrated
!>        with (for a displaced Lagrangian parcel, the zone of the cell it comes from). To be used with
!>        a grouping by OLD zones (set_wat_indices_zones on solid_chemistry_old). Previous-step
!>        reaction term: kinetic only (the equilibrium block cancels under U_z).
! subroutine mix_within_zone_old_rz(all_waters,zone_mix_idx,ratios_conc_norm,ratios_reacts_norm,&
!     Delta_t,theta_r,zone_water,ideal,rep_rz,zone_solid_out)
!     implicit none
!     type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
!     integer(kind=4), intent(in) :: zone_mix_idx(:) !> global indices of the group's mixing waters (all with OLD zone z)
!     real(kind=8), intent(in) :: ratios_conc_norm(:) !> within-zone concentration mixing ratios (sum=1), aligned with zone_mix_idx
!     real(kind=8), intent(in) :: ratios_reacts_norm(:) !> within-zone reaction mixing ratios (sum=1), aligned with zone_mix_idx
!     real(kind=8), intent(in) :: Delta_t !> time step
!     real(kind=8), intent(in) :: theta_r !> time weighting factor
!     type(aqueous_chemistry_c), intent(out) :: zone_water !> temporary water of zone z holding the speciated sub-mixture
!     logical, intent(in) :: ideal !> .true. ideal speciation; .false. non-ideal
!     class(reactive_zone_c), intent(in), target :: rep_rz !> reactive zone defining zone z (representative's OLD or CURRENT reactive zone); target: the scratch solid chemistry points to it
!     type(solid_chemistry_c), intent(out) :: zone_solid_out !> copy of the sub-mixture's scratch solid chemistry

!     call mix_within_zone_core(all_waters,zone_mix_idx,ratios_conc_norm,ratios_reacts_norm,Delta_t,theta_r,&
!         zone_water,ideal,rep_rz,zone_solid_out)
! end subroutine mix_within_zone_old_rz

!> \brief Within-zone mixing in the NEW (current) reactive zone of the group
!>        (solid_chemistry%reactive_zone of its representative): the zone of the cell the members sit
!>        in now. To be used with a grouping by CURRENT zones (set_wat_indices_zones on
!>        solid_chemistry): the displaced parcels are brought to equilibrium with their new zone by
!>        this speciation. Previous-step reaction term: kinetic only (the equilibrium block cancels
!>        under U_z).
subroutine mix_within_zone_new_rz(this,all_waters,zone_mix_idx,ratios_conc_norm,ratios_reacts_norm,&
    Delta_t,theta_r,ideal,rep_rz,rep_mz)
    implicit none
    class(aqueous_chemistry_c), intent(inout) :: this !> temporary water of zone z holding the speciated sub-mixture
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: zone_mix_idx(:) !> global indices of the group's mixing waters (all with CURRENT zone z)
    real(kind=8), intent(in) :: ratios_conc_norm(:) !> within-zone concentration mixing ratios (sum=1), aligned with zone_mix_idx
    real(kind=8), intent(in) :: ratios_reacts_norm(:) !> within-zone reaction mixing ratios (sum=1), aligned with zone_mix_idx
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta_r !> time weighting factor
    logical, intent(in) :: ideal !> .true. ideal speciation; .false. non-ideal
    class(reactive_zone_c), intent(in), target :: rep_rz !> reactive zone defining zone z (representative's OLD or CURRENT reactive zone); target: the scratch solid chemistry points to it
    class(mineral_zone_c), intent(in), target :: rep_mz !> mineral zone defining zone z (representative's OLD or CURRENT reactive zone); target: the scratch solid chemistry points to it

    call mix_within_zone_core(this,all_waters,zone_mix_idx,ratios_conc_norm,ratios_reacts_norm,Delta_t,theta_r,&
        ideal,rep_rz,rep_mz)
end subroutine mix_within_zone_new_rz

!> \brief Consistent WMA cross-zone mixing (algorithm stages 2 and 3) for one target cell.
!> \details Stage 2: mixes, within each reactive zone, the target waters of that zone's column of
!>          wat_indices_zones (the target waters grouped by the zone of their previous-step state:
!>          the zone their old concentrations were equilibrated with, i.e. for a displaced Lagrangian
!>          parcel the zone of the cell it comes from) with their mixing ratios in this cell's stencil
!>          (mix_within_zone) into a temporary
!>          speciated water, each water carrying its reaction mixing ratio in the within-zone kinetic
!>          term. Stage 3: cross-zone mixes those zone waters in AQUEOUS-PHASE species space
!>          (so different primary/secondary partitions across zones are handled through each zone
!>          water's indices_aq_species), weighted by each zone's total concentration mixing ratio,
!>          together with the external waters (ext_wat_indices, mixed directly):
!>            c_hat = sum_z lambda_c,z c_hat_z + sum_ext lambda_ext c_ext,
!>          while c_mix is the purely conservative mix of all the waters (no reaction).
!>          Direct mix (no internal zone in the stencil, single-zone domain, or every non-external water
!>          from the target's own current zone): mix_within_zone is skipped and the (1-theta_r) previous-step
!>          KINETIC term of the reaction-mixing waters is added here to c_hat (c_mix stays conservative).
!>          The new-step (theta_r) reaction term, c_hat + Delta_t*r_tilde_i, is added afterwards by the
!>          solver (add_curr_step_react_mix_contrib_opt* already applies the theta_r weight to r_tilde_i);
!>          the self kinetic new-step term is added by the Newton solver.
subroutine compute_c_mix_by_zones_cons(all_waters,zone_waters,conc_idx,conc_lam,react_lam,&
    wat_indices_zones,ext_wat_indices,Delta_t,theta_r,c_hat,c_mix,ideal)
    implicit none
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:) !> full waters array (self water's conc_mix is updated with the conservative mix)
    type(zone_water_c), intent(inout) :: zone_waters(:) !> one temporary speciated water per zone
    integer(kind=4), intent(in) :: conc_idx(:) !> concentration-mixing waters (self first), global indices
    real(kind=8), intent(in) :: conc_lam(:) !> concentration mixing ratios (aligned with conc_idx, sum=1)
    real(kind=8), intent(in) :: react_lam(:) !> reaction mixing ratios (aligned with conc_idx; 0 for the waters that mix concentrations but not reaction rates)
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (col j = zone id): global indices of the target waters whose previous-step state belongs to that zone
    integer(kind=4), intent(in) :: ext_wat_indices(:) !> global indices of the external (boundary and recharge) waters
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta_r !> time weighting factor
    real(kind=8), intent(out), allocatable :: c_hat(:) !> cross-zone mixed AND within-zone reacted concentrations (target basis); allocated here, size = num_var_act_species
    real(kind=8), intent(out), allocatable :: c_mix(:) !> conservative part of c_hat (aqueous variable-activity species only); allocated here, size = num_aq_var_act_species
    logical, intent(in) :: ideal !> per-zone speciation: .true. ideal (gamma=1), .false. non-ideal

    integer(kind=4) :: p, n_vas
    logical :: direct_mix
    real(kind=8), allocatable :: c_hat_zones(:), c_mix_zones(:) !> stages 2+3 result (cross-zone mix of the speciated zone sub-mixtures), before the kinetic term
    real(kind=8), allocatable :: Sk_nc_old_j(:,:), rk_old_j(:), kin_prev(:)

    call cross_zone_mix_core(all_waters,zone_waters,conc_idx,conc_lam,react_lam,wat_indices_zones,ext_wat_indices,&
        Delta_t,theta_r,ideal,c_hat_zones,c_mix_zones,direct_mix)
    c_hat=c_hat_zones
    c_mix=c_hat
    !> Direct mix: the within-zone kinetic term was not formed, so add the reaction-mixing waters'
    !> (1-theta_r) previous-step kinetic term here; in the multi-zone case mix_within_zone carries it.
    if (direct_mix) then
        n_vas=size(c_hat_zones)
        allocate(kin_prev(n_vas))
        kin_prev=0d0
        do p=1,size(conc_idx)
            if (react_lam(p)==0d0) cycle
            if (all_waters(conc_idx(p))%indices_rk_old%num_cols>0) then
                Sk_nc_old_j=all_waters(conc_idx(p))%get_Sk_nc_old()
                rk_old_j=all_waters(conc_idx(p))%get_rk_old()
                if (size(Sk_nc_old_j,2)==n_vas) &
                    kin_prev=kin_prev+react_lam(p)*matmul(transpose(Sk_nc_old_j),rk_old_j)
            end if
        end do
        call add_limited_increment(c_mix,Delta_t*(1d0-theta_r)*kin_prev,c_hat)
    end if
    call write_zone_dbg_c_mix(all_waters(conc_idx(1)),c_hat,c_mix)
end subroutine compute_c_mix_by_zones_cons

!> \brief Lumped WMA cross-zone mixing for one target cell.
!> \details Same stages 2 and 3 as compute_c_mix_by_zones_cons but the within-zone mixes carry NO
!>          kinetic term (equilibrium-only sub-mixtures, doc eq. 2); the previous-step kinetic
!>          contribution of the target water alone is added after the cross-zone mix, weighted by the
!>          lumped reaction mixing ratio (doc eq. 4):  c_hat += lambda_L Delta_t (1-theta_r) Sk_i^T rk_i^k.
subroutine compute_c_mix_by_zones_lump(all_waters,zone_waters,conc_idx,conc_lam,lambda_L,&
    wat_indices_zones,ext_wat_indices,Delta_t,theta_r,c_hat,c_mix,ideal)
    implicit none
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:) !> full waters array (self water's conc_mix is updated with the conservative mix)
    type(zone_water_c), intent(inout) :: zone_waters(:) !> one temporary speciated water per zone
    integer(kind=4), intent(in) :: conc_idx(:) !> concentration-mixing waters (self first), global indices
    real(kind=8), intent(in) :: conc_lam(:) !> concentration mixing ratios (aligned with conc_idx, sum=1)
    real(kind=8), intent(in) :: lambda_L !> lumped reaction mixing ratio of the target water
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (col j = zone id): global indices of the waters whose previous-step state belongs to that zone
    integer(kind=4), intent(in) :: ext_wat_indices(:) !> global indices of the external (boundary and recharge) waters
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta_r !> time weighting factor
    real(kind=8), intent(out), allocatable :: c_hat(:) !> cross-zone mixed AND within-zone reacted concentrations (target basis); allocated here, size = num_var_act_species
    real(kind=8), intent(out), allocatable :: c_mix(:) !> conservative part of c_hat (aqueous variable-activity species only); allocated here, size = num_aq_var_act_species
    logical, intent(in) :: ideal !> per-zone speciation: .true. ideal (gamma=1), .false. non-ideal

    logical :: direct_mix
    real(kind=8), allocatable :: react_lam(:) !> within-zone reaction ratios: all zero (equilibrium-only sub-mixtures)
    real(kind=8), allocatable :: c_hat_zones(:) !> stages 2+3 result: cross-zone mix of the speciated zone sub-mixtures (full var-act layout), before the kinetic term
    real(kind=8), allocatable :: c_mix_zones(:) !> conservative aqueous part of c_hat_zones (pure previous-step mix)
    real(kind=8), allocatable :: Sk_nc_old(:,:), rk_old(:)
    logical, save :: mix_dbg_done=.false. !> [MIX-DBG] one-shot flag of this call site

    allocate(react_lam(size(conc_idx)))
    react_lam=0d0
    if (within_zone_dbg_step==1 .and. .not. mix_dbg_done) then
        mix_dbg_done=.true.
        call print_mix_dbg_stencil("compute_c_mix_by_zones_lump",conc_idx,conc_lam,lambda_L)
    end if
    call cross_zone_mix_core(all_waters,zone_waters,conc_idx,conc_lam,react_lam,wat_indices_zones,ext_wat_indices,&
        Delta_t,theta_r,ideal,c_hat_zones,c_mix_zones,direct_mix)
    !> Both outputs are set unconditionally: c_hat must be allocated also when the target has no
    !> kinetic reactions (equilibrium-only runs), since the solver uses it right away
    c_hat=c_hat_zones
    c_mix=c_hat
    associate(self => all_waters(conc_idx(1)))
    if (self%indices_rk_old%num_cols>0) then
        Sk_nc_old=self%get_Sk_nc_old()
        rk_old=self%get_rk_old()
        if (size(Sk_nc_old,2)==size(c_hat_zones)) &
            call add_limited_increment(c_mix,Delta_t*(1d0-theta_r)*lambda_L*matmul(transpose(Sk_nc_old),rk_old),c_hat)
    end if
    call write_zone_dbg_c_mix(self,c_hat,c_mix)
    end associate
end subroutine compute_c_mix_by_zones_lump

!> [MIX-DBG] stdout trace of one target cell's mixing stencil: concentration-mixing waters (self
!> first), their concentration mixing ratios (column-aligned with the indices), optionally the
!> reaction mixing ratios (same alignment) and the lumped reaction mixing ratio, and the
!> concentration ratio sum (must be 1). Each call site gates it on the first time step
!> (within_zone_dbg_step==1) and prints only once (first target water processed in that step).
subroutine print_mix_dbg_stencil(mixer,conc_idx,conc_lam,lambda_L,react_lam)
    implicit none
    character(len=*), intent(in) :: mixer !> name of the calling mixer
    integer(kind=4), intent(in) :: conc_idx(:) !> concentration-mixing waters (self first), global indices
    real(kind=8), intent(in) :: conc_lam(:) !> concentration mixing ratios aligned with conc_idx
    real(kind=8), intent(in), optional :: lambda_L !> lumped reaction mixing ratio of the target water
    real(kind=8), intent(in), optional :: react_lam(:) !> reaction mixing ratios aligned with conc_idx
    if (present(lambda_L)) then
        print '(A,I0,A,I0,A,I0,A,ES12.4,A,A)', "[MIX-DBG] step ", within_zone_dbg_step, " target water=", conc_idx(1), &
            " nmix=", size(conc_idx), " lambda_L=", lambda_L, " mixer=", mixer
    else
        print '(A,I0,A,I0,A,I0,A,A)', "[MIX-DBG] step ", within_zone_dbg_step, " target water=", conc_idx(1), &
            " nmix=", size(conc_idx), " mixer=", mixer
    end if
    print '(("[MIX-DBG]   conc_idx  =",12(1X,I12)))', conc_idx
    print '(("[MIX-DBG]   conc_lam  =",12(1X,ES12.4)))', conc_lam
    if (present(react_lam)) print '(("[MIX-DBG]   react_lam =",12(1X,ES12.4)))', react_lam
    print '(A,ES23.15)', "[MIX-DBG]   sum(conc_lam) =", sum(conc_lam)
end subroutine print_mix_dbg_stencil

!> Patankar-type scalar limiter: c_hat = c_mix + omega*delta with a uniform omega in [0,1] so that no
!> species is driven below floor_frac of its value; a uniform omega keeps omega*Sk^T*rk a valid reaction extent.
subroutine add_limited_increment(c_mix,delta,c_hat)
    implicit none
    real(kind=8), intent(in) :: c_mix(:)
    real(kind=8), intent(in) :: delta(:)
    real(kind=8), intent(out) :: c_hat(:)
    integer(kind=4) :: iv
    real(kind=8) :: omega
    real(kind=8), parameter :: floor_frac=1d-3 !> keep at least this fraction of each species
    omega=1d0
    do iv=1,size(c_mix)
        if (delta(iv)<0d0 .and. c_mix(iv)>0d0) &
            omega=min(omega,(1d0-floor_frac)*c_mix(iv)/(-delta(iv)))
    end do
    if (omega<0d0) omega=0d0
    c_hat=c_mix+omega*delta
end subroutine add_limited_increment

!> [ZONE-DBG] final rows of the traced cell: c_hat (aqueous primaries then aqueous secondaries, skipping
!> the free site), the conservative mix c_mix and their difference; then disarm the trace.
subroutine write_zone_dbg_c_mix(self,c_hat,c_mix)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: self
    real(kind=8), intent(in) :: c_hat(:), c_mix(:)
    integer(kind=4) :: jj, n_aq_prim, n_aq_sec, n_prim
    real(kind=8), allocatable :: chat_aq(:)
    if (within_zone_dbg_active) then
        n_aq_prim=size(self%ind_prim_species)
        n_aq_sec=self%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
        n_prim=self%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        allocate(chat_aq(n_aq_prim+n_aq_sec))
        chat_aq(1:n_aq_prim)=c_hat(1:n_aq_prim)
        chat_aq(n_aq_prim+1:n_aq_prim+n_aq_sec)=c_hat(n_prim+1:n_prim+n_aq_sec)
        write(within_zone_dbg_unit,'("[ZONE-DBG]   aq var act species =",*(A13))') &
            (trim(self%aq_phase%aq_species(self%indices_aq_phase(self%ind_var_act_species(jj)))%name), jj=1,n_aq_prim), &
            (trim(self%aq_phase%aq_species(self%indices_aq_phase(self%ind_var_act_species(n_prim+jj)))%name), jj=1,n_aq_sec)
        write(within_zone_dbg_unit,'("[ZONE-DBG]   c_hat              =",*(ES13.5))') chat_aq
        write(within_zone_dbg_unit,'("[ZONE-DBG]   c_mix              =",*(ES13.5))') c_mix
        write(within_zone_dbg_unit,'("[ZONE-DBG]   within-zone rxn    =",*(ES13.5))') chat_aq-c_mix
    end if
    within_zone_dbg_active=.false. !> [ZONE-DBG] disarm the trace for the next (possibly untraced) cell
end subroutine write_zone_dbg_c_mix

!> Conservative core shared by compute_c_mix_by_zones_{cons,lump}: stage 2 (within-zone mixing with
!> the per-water reaction ratios react_lam) and stage 3 (cross-zone aqueous mixing). The grouping
!> of the mixing waters is taken as given, no classification is done here: the members of reactive
!> zone z are the target waters listed in column z of wat_indices_zones (the zone their
!> previous-step state was equilibrated with: for a displaced Lagrangian parcel the zone of the cell
!> it comes from, the solvers regroup the columns after each displacement; with a Neumann outflow the
!> boundary water holding the parcel that left the last cell is listed in that cell's zone column too)
!> and the external waters
!> are those of ext_wat_indices; each one's mixing ratio is read by its global water index. The zone
!> members are mixed and re-speciated within their zone (mix_within_zone), whose component matrix
!> cancels their equilibrium reaction term (U_z*S_eq^T = 0); the external waters (boundary and
!> recharge, whose water-type zones are not in chemistry%reactive_zones) are mixed directly in
!> stage 3; the cross-zone mix is then equilibrated with the target cell's CURRENT zone by the solver.
!> Builds c_hat (mixed + within-zone reacted) and c_mix (conservative mix) and stores the latter on
!> the target water. Adds NO kinetic term itself: the caller adds its own previous-step kinetic
!> contribution (direct_mix tells it whether mix_within_zone was skipped).
subroutine cross_zone_mix_core(all_waters,zone_waters,conc_idx,conc_lam,react_lam,wat_indices_zones,&
    ext_wat_indices,Delta_t,theta_r,ideal,c_hat,c_mix,direct_mix)
    implicit none
    type(aqueous_chemistry_c), intent(inout) :: all_waters(:) !> full waters array (self water's conc_mix is updated with the conservative mix)
    type(zone_water_c), intent(inout) :: zone_waters(:) !> one temporary speciated water per zone
    integer(kind=4), intent(in) :: conc_idx(:) !> concentration-mixing waters (self first), global indices
    real(kind=8), intent(in) :: conc_lam(:) !> concentration mixing ratios (aligned with conc_idx, sum=1)
    real(kind=8), intent(in) :: react_lam(:) !> reaction mixing ratio of each concentration-mixing water (aligned with conc_idx; 0 for the waters that mix concentrations but not reaction rates) for the within-zone kinetic term
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (col j = zone id): global indices of the waters whose previous-step state belongs to that zone (chemistry%wat_indices_zones): the members of each within-zone mix
    integer(kind=4), intent(in) :: ext_wat_indices(:) !> global indices of the external waters (chemistry%ext_waters_indices: boundary and recharge; the Lagrangian solvers add the inflow cells, whose parcel is the boundary water, and remove the Neumann outflow waters, which hold the parcel that left the last cell and are grouped in that cell's zone column instead): mixed directly in stage 3, never speciated within a zone
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta_r !> time weighting factor
    logical, intent(in) :: ideal !> per-zone speciation: .true. ideal (gamma=1), .false. non-ideal
    real(kind=8), intent(out), allocatable :: c_hat(:) !> cross-zone mixed AND within-zone reacted concentrations (target basis); allocated here, size = num_var_act_species (aqueous + non-aqueous variable-activity species)
    real(kind=8), intent(out), allocatable :: c_mix(:) !> conservative part of c_hat: pure previous-step mixing, no reaction rates; allocated here, size = num_aq_var_act_species (aqueous variable-activity species only)
    logical, intent(out) :: direct_mix !> .true. when mix_within_zone was skipped (all waters mixed directly): no contributing zone other than the target's own current zone, or a single-zone domain

    integer(kind=4) :: nmix, p, q, w, z, nzd, self_zone_id, i, aq_idx
    integer(kind=4) :: n_aq, n_aq_prim, n_aq_sec, n_prim, n_vas
    real(kind=8), allocatable :: lam_glob(:) !> concentration mixing ratio of every water, by global index (0 outside this cell's stencil)
    real(kind=8), allocatable :: rlam_glob(:) !> reaction mixing ratio of every water, by global index (0 outside this cell's stencil)
    real(kind=8), allocatable :: lam_c_zone(:) !> per-zone total concentration mixing ratio, indexed by zone id (0 for a zone absent from the stencil)
    real(kind=8) :: lam_ext !> total concentration mixing ratio of the external waters
    real(kind=8), allocatable :: grp_rc(:), grp_rr(:) !> group within-zone conc & reaction ratios
    real(kind=8) :: sum_rr
    
    !type(aqueous_chemistry_c), allocatable :: zone_waters_aq_chem(:) !> aqueous chemistry of the zone waters
    real(kind=8), allocatable :: c_tmp(:)
    real(kind=8), allocatable :: c_cons(:) !> cross-zone aqueous mix before speciation, aq_phase order: direct mix = pure previous-step mix; multi-zone = weighted sum of the zone waters' pre-speciation mixes c_hat_z (conc_mix), i.e. including the limited previous-step kinetic term
    logical, parameter :: dbg=.false. !> set .true. to enable per-reactive-zone mixing debug prints
    logical, save :: zonedbg_descr_written = .false. !> [ZONE-DBG] one-time file legend written flag
    logical, save :: mix_dbg_done=.false. !> [MIX-DBG] one-shot flag of this call site
    logical :: single_zone !> .true. when the whole domain has a single reactive zone (per-zone lambda renormalization is then skipped)

    nmix=size(conc_idx) !> number of mixing waters
    nzd=wat_indices_zones%num_cols !> number of domain reactive zones (zone id = column index)
    single_zone=(nzd==1) !> flag indicating if the whole domain has a single reactive zone
    !> [MIX-DBG] stencil of the first target water of the first time step: concentration and
    !> reaction mixing ratios of the mixing waters, as received by the core
    if (within_zone_dbg_step==1 .and. .not. mix_dbg_done) then
        mix_dbg_done=.true.
        call print_mix_dbg_stencil("cross_zone_mix_core",conc_idx,conc_lam,react_lam=react_lam)
    end if
    !> Mixing ratios by GLOBAL water index (the stencil scattered over all_waters, 0 elsewhere), so the
    !> zone groups and the external waters are read straight from wat_indices_zones and
    !> ext_wat_indices
    allocate(lam_glob(size(all_waters)),rlam_glob(size(all_waters)))
    lam_glob=0d0
    rlam_glob=0d0
    lam_glob(conc_idx)=conc_lam
    rlam_glob(conc_idx)=react_lam
    allocate(lam_c_zone(nzd))
    do z=1,nzd
        lam_c_zone(z)=sum(lam_glob(wat_indices_zones%cols(z)%col_1))
    end do
    lam_ext=sum(lam_glob(ext_wat_indices))
    !> Direct mix (mix_within_zone skipped, all waters mixed directly; the caller adds the kinetic
    !> term) when no zone other than the TARGET's own CURRENT reactive zone contributes to the
    !> stencil (only external waters and/or waters from the target's zone), or when the domain has
    !> a single reactive zone. In that case the target's equilibration conserves exactly the
    !> components the within-zone speciation would (U_z*c_hat_z = U_z*c_mix_z), so the stage-2
    !> speciation is redundant. A contributing OLD zone DIFFERENT from the target's current zone
    !> (parcels displaced across a zone interface together with their stencil) still needs the
    !> within-zone speciation in that old zone.
    !> the target's CURRENT reactive zone (the within-zone mixing uses the new-zone variant)
    self_zone_id=all_waters(conc_idx(1))%solid_chemistry%reactive_zone%id
    direct_mix=.true.
    do z=1,nzd
        if (z/=self_zone_id .and. lam_c_zone(z)>0d0) direct_mix=.false.
    end do
    if (single_zone) direct_mix=.true.

    !> [ZONE-DBG] arm the per-cell zone-mixing trace for this target water. mix_within_zone reads
    !> within_zone_dbg_active to decide whether to emit its per-zone rows; disarmed at routine exit.
    !> Only the selected cells (within_zone_dbg_cells) during the first within_zone_dbg_max_step
    !> steps are traced, keeping the file small (like the CELL-DBG trace).
    within_zone_dbg_active=.false.
    if (within_zone_dbg .and. within_zone_dbg_unit/=-1 .and. within_zone_dbg_step>=1 .and. &
        within_zone_dbg_step<=within_zone_dbg_max_step .and. allocated(within_zone_dbg_cells)) then
        within_zone_dbg_active=any(within_zone_dbg_cells==conc_idx(1))
    end if
    if (within_zone_dbg_active .and. .not. zonedbg_descr_written) then
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] =================================================================="
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] REMIX per-cell zone-mixing trace (compute_c_mix_by_zones / mix_within_zone)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] one block per (step, target water); concentrations [mol/L]"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] stage 2 (mix_within_zone) emits one sub-block per reactive zone in the"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   cell's mixing stencil:"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   zone var act : the zone's variable-activity species (conc_nc column order)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   ratios_conc  : within-zone concentration mixing ratios (per member, sum=1)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   r_tilde_z    : within-zone previous-step kinetic reaction term"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   c_mix_z      : conservative within-zone mix entering speciation"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   c_hat_z      : speciated within-zone sub-mixture (Newton output)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   rho_eq/Re/re_mean : eq contributions, amounts and rates (aqueous & solid)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] stage 3 (cross-zone aqueous mixing), target-cell aqueous var-act species:"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   lambda_c per zone : each zone's total concentration mixing ratio"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   c_hat        : cross-zone mixed (and within-zone reacted) var-act conc = result of this subroutine"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   c_mix        : conservative cross-zone mix (transport only, no reaction)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG]   within-zone rxn : c_hat - c_mix (per aqueous var-act species)"
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] =================================================================="
        zonedbg_descr_written=.true.
    end if
    if (within_zone_dbg_active) then
        write(within_zone_dbg_unit,'(A)') "[ZONE-DBG] ------------------------------------------------------------------"
        write(within_zone_dbg_unit,'(A,I0,A,I0,A,I0,A,I0,A,L1)') "[ZONE-DBG] step=", within_zone_dbg_step, &
            " target water index=", conc_idx(1), " nmix=", nmix, " nzones=", count(lam_c_zone>0d0), " direct_mix=", direct_mix
        write(within_zone_dbg_unit,'("[ZONE-DBG]   zone ids =",*(1X,I0))') pack([(z,z=1,nzd)],lam_c_zone>0d0)
    end if

    if (dbg) then
        print '(A,I0,A,I0,A,I0)', "[ZONE-MIX] self=", conc_idx(1), " nmix=", nmix, " nzones=", count(lam_c_zone>0d0)
        print '(A,*(I0,1X))', "[ZONE-MIX]   zone_ids=", pack([(z,z=1,nzd)],lam_c_zone>0d0)
        print '(A,*(I0,1X))', "[ZONE-MIX]   mix_glob=", conc_idx
        print '(A,*(ES12.4,1X))', "[ZONE-MIX]   lam_zone=", lam_c_zone
        print '(A,ES12.4)', "[ZONE-MIX]   lam_ext =", lam_ext
    end if

    !> Stage 2: within-zone mixing -> one temporary speciated water per zone (multi-zone domains only),
    !> indexed by zone id. The members of zone z are the waters of column z of
    !> wat_indices_zones, each with its own mixing ratios (0 for a member outside the stencil). A
    !> zone absent from the stencil, or present with ZERO total weight (a far-field cell of the dense
    !> implicit theta_t=1 stencil), contributes nothing: it is skipped here (avoiding a 0/0 in the
    !> ratio renormalisation) and in stage 3 via the same lam_c_zone(z)<=0 test.
    !> Every zone water of this target starts the step as "not built" (lambda_c=0); only the zones
    !> built below get lambda_c>0, which is how write_zone_waters (chemistry_m) picks the current
    !> step's sub-mixtures out of chemistry%zone_waters.
    do z=1,size(zone_waters)
        call zone_waters(z)%set_lambda_c(0d0)
    end do
    !allocate(zone_waters(nzd))
    if (.not. direct_mix) then
    do z=1,nzd
        if (lam_c_zone(z)<=0d0) cycle
        associate(grp_glob => wat_indices_zones%cols(z)%col_1)
        grp_rc=lam_glob(grp_glob)
        grp_rr=rlam_glob(grp_glob)
        !> Per-zone lambda renormalization applies only with MULTIPLE reactive zones. With a single
        !> reactive zone the raw mixing ratios are used directly (no rescaling), even when external
        !> waters are present; the speciated sub-mixture then carries its own weight (stage 3 uses
        !> weight 1 for it).
        !if (.not. single_zone) then
            grp_rc=grp_rc/lam_c_zone(z) !> within-zone concentration ratios (sum=1)
            sum_rr=sum(grp_rr)
            if (sum_rr>0d0) grp_rr=grp_rr/sum_rr !> within-zone reaction ratios (sum=1)
        !end if
        if (dbg) print '(A,I0,A,I0,A,ES12.4)', "[ZONE-MIX]   zone_id=", z, &
            " n_members=", size(grp_glob), " lambda_c_zone=", lam_c_zone(z)
        if (within_zone_dbg_active) write(within_zone_dbg_unit,'(A,I0,A,ES13.5)') &
            "[ZONE-DBG]   -- reactive zone id=", z, " lambda_c=", lam_c_zone(z)
        call zone_waters(z)%aq_chem%mix_within_zone_new_rz(all_waters,grp_glob,grp_rc,grp_rr,Delta_t,theta_r,ideal,&
            all_waters(conc_idx(1))%solid_chemistry%reactive_zone,all_waters(conc_idx(1))%solid_chemistry%mineral_zone) !> rep_rz = self water's zone
        call zone_waters(z)%set_tar_wat_idx(conc_idx(1))
        call zone_waters(z)%set_zone_id(z)
        call zone_waters(z)%set_lambda_c(lam_c_zone(z))
        end associate
    end do
    end if

    !> Stage 3: cross-zone conservative mixing in aqueous-phase species space.
    associate(self => all_waters(conc_idx(1)))
    n_aq      = self%aq_phase%num_species
    n_aq_prim = size(self%ind_prim_species)
    n_aq_sec  = self%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
    n_prim    = self%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_vas     = self%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(c_tmp(n_aq)); c_tmp=0d0
    allocate(c_cons(n_aq)); c_cons=0d0
    do i=1,n_aq
        aq_idx=self%indices_aq_phase(i)
        if (direct_mix) then
            !> Direct mix: mix_within_zone was skipped (zone_waters not built), so the aqueous slots are
            !> the purely conservative previous-step mix over ALL mixing waters (domain and external,
            !> each mapped through its own aqueous-species index); the kinetic previous-step term is
            !> added by the caller.
            do p=1,nmix
                c_tmp(i)=c_tmp(i)+conc_lam(p)*all_waters(conc_idx(p))%conc_old(&
                    all_waters(conc_idx(p))%indices_aq_species(aq_idx))
            end do
            c_cons(i)=c_tmp(i)
            cycle
        end if
        !> domain zones: speciated within-zone sub-mixtures. c_tmp takes the speciated concentration;
        !> c_cons the pre-speciation mix c_hat_z, read from the zone water's conc_mix through the
        !> species' slot (var-act column order, NOT the water's concentration index: see
        !> zone_water_mix_cons).
        do z=1,nzd
            if (lam_c_zone(z)<=0d0) cycle !> zone absent from the stencil or with zero weight: zone_waters(z) not built
            if (single_zone) then
                !> single reactive zone: raw (un-renormalized) ratios were used, so the sub-mixture
                !> already carries its own weight -> combine with weight 1 (not lam_c_zone).
                c_tmp(i)=c_tmp(i)+zone_waters(z)%aq_chem%concentrations(zone_waters(z)%aq_chem%indices_aq_species(aq_idx))
                c_cons(i)=c_cons(i)+zone_water_mix_cons(zone_waters(z)%aq_chem,zone_waters(z)%aq_chem%indices_aq_species(aq_idx))
            else
                c_tmp(i)=c_tmp(i)+lam_c_zone(z)*zone_waters(z)%aq_chem%concentrations(&
                    zone_waters(z)%aq_chem%indices_aq_species(aq_idx))
                c_cons(i)=c_cons(i)+lam_c_zone(z)*zone_water_mix_cons(zone_waters(z)%aq_chem,&
                    zone_waters(z)%aq_chem%indices_aq_species(aq_idx))
            end if
        end do
        !> external waters: previous-step aqueous concentration mixed directly (mapped through each
        !> water's own aqueous-species index), not per-zone speciated
        do q=1,size(ext_wat_indices)
            w=ext_wat_indices(q)
            c_tmp(i)=c_tmp(i)+lam_glob(w)*all_waters(w)%conc_old(all_waters(w)%indices_aq_species(aq_idx))
            c_cons(i)=c_cons(i)+lam_glob(w)*all_waters(w)%conc_old(all_waters(w)%indices_aq_species(aq_idx))
        end do
    end do
    !> Mass-balance check. (1) Partition: the mixing weights must sum to 1 (conservative transport).
    !> (2) Per-species (dbg): compare the mixed aqueous concentration against the purely conservative
    !> cross-zone mix of the previous-step concentrations; the difference is the within-zone reaction
    !> (~0 for inert species such as Cl-, nonzero for reacting species).
    block
        real(kind=8) :: wsum, wsum_raw
        integer(kind=4) :: zz
        !> Invariant: the concentration mixing ratios of this mixing zone (the target cell's mixing
        !> stencil) must sum to 1 (conservative transport). Check the raw ratios directly, and
        !> cross-check the per-reactive-zone partition (sum of per-zone totals + external total) so
        !> a stencil water listed neither in a zone column nor as external is caught too.
        wsum_raw=sum(conc_lam(1:nmix))
        if (direct_mix) then
            wsum=wsum_raw !> direct mix: all waters mixed conservatively (weights sum to 1)
        else
            wsum=sum(lam_c_zone)+lam_ext
        end if
        if (abs(wsum_raw-1d0)>1d-10 .or. abs(wsum-1d0)>1d-10) then
            print '(A,I0,A,ES13.5,A,ES13.5)', &
                "[ZONE-MIX MB] ERROR: self=", conc_idx(1), " sum(mixing ratios)=", wsum_raw, &
                " per-zone partition sum=", wsum
            do zz=1,nzd
                print '(A,I0,A,ES13.5)', "[ZONE-MIX MB]   reactive zone id=", zz, &
                    " sum(mixing ratios)=", lam_c_zone(zz)
            end do
            print '(A,ES13.5)', "[ZONE-MIX MB]   external waters sum(mixing ratios)=", lam_ext
            error stop "[ZONE-MIX MB] mixing ratios do not sum to 1 in the mixing zone"
        end if
    end block
    !> c_hat: cross-zone mixed (and within-zone reacted) concentration over ALL variable-activity
    !> species (allocated here, size num_var_act_species). Aqueous slots come from the zone waters; the
    !> non-aqueous variable-activity tail (solid/gas primary, surface complexes, gas eq var act) is
    !> taken from the target water's previous step via get_conc_old_var_act_new_rz. Full layout [aq prim]
    !> [solid prim][gas prim][aq sec var act][surf compl][gas eq var act], so aqueous secondaries sit
    !> at offset n_prim.
    if (allocated(c_hat)) deallocate(c_hat)
    allocate(c_hat(n_vas))
    c_hat=self%get_conc_old_var_act_new_rz()
    do i=1,n_aq_prim
        c_hat(i)=c_tmp(self%ind_prim_species(i))
    end do
    do i=1,n_aq_sec
        c_hat(n_prim+i)=c_tmp(self%ind_sec_species(i))
    end do
    if (dbg) print '(A,I0,A,*(ES12.4))', "[ZONE-MIX] self=", conc_idx(1), " c_hat=", c_hat
    !> c_mix: conservative part of c_hat (allocated here, size num_aq_var_act_species): AQUEOUS
    !> variable-activity species ONLY, compact layout [aqueous primary][aqueous secondary var act]
    !> (aqueous secondaries at offset n_aq_prim), from the pure previous-step mix.
    if (allocated(c_mix)) deallocate(c_mix)
    allocate(c_mix(n_aq_prim+n_aq_sec))
    c_mix=0d0
    do i=1,n_aq_prim
        c_mix(i)=c_cons(self%ind_prim_species(i))
    end do
    do i=1,n_aq_sec
        c_mix(n_aq_prim+i)=c_cons(self%ind_sec_species(i))
    end do
    !> Store the conservative aqueous mix on the self (target) water for the conc_mix output variable
    !> (compact [aq prim | aq sec] layout, as compute_c_mix_cons_global does); the speciation input
    !> c_hat (full layout, plus the new-step reaction term) is stored by the solver in conc_hat.
    self%conc_mix=c_mix
    if (within_zone_dbg_active .and. .not. direct_mix) &
        write(within_zone_dbg_unit,'("[ZONE-DBG]   lambda_c per zone  =",*(ES13.5))') lam_c_zone
    end associate
end subroutine cross_zone_mix_core

!> Pre-speciation mix value of a zone water for one of its aqueous species. conc_hat is laid out
!> [aq prim (1:n_p) | aq sec var act (n_p+1:n_p+n_s)] in the zone's var-act COLUMN order, while k is
!> the species' index in the water's concentrations array, so the slot is found through
!> ind_prim_species / ind_sec_species (which hold those concentration indices in column order, see
!> set_ind_species_react_zone). A species without a slot (constant-activity aqueous species such as
!> water, which is neither mixed nor reacted) returns its speciated concentration.
function zone_water_mix_cons(zw,k) result(val)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: zw !> zone water (speciated sub-mixture)
    integer(kind=4), intent(in) :: k !> index of the species in zw%concentrations
    real(kind=8) :: val
    integer(kind=4) :: j,n_p,n_s
    val=zw%concentrations(k)
    if (.not. allocated(zw%conc_hat) .or. .not. allocated(zw%ind_prim_species)) return
    n_p=min(size(zw%ind_prim_species),size(zw%conc_hat))
    do j=1,n_p
        if (zw%ind_prim_species(j)==k) then
            val=zw%conc_hat(j)
            return
        end if
    end do
    if (.not. allocated(zw%ind_sec_species)) return
    n_s=min(size(zw%ind_sec_species),size(zw%conc_hat)-n_p)
    do j=1,n_s
        if (zw%ind_sec_species(j)==k) then
            val=zw%conc_hat(n_p+j)
            return
        end if
    end do
end function zone_water_mix_cons

!> Appends one within-zone sub-mixture to the current-step output store (see zone_water_out_c).
subroutine set_zone_water_attributes(this,tar_wat_idx,zone_id,lambda_c,aq_chem)
    implicit none
    class(zone_water_c), intent(inout) :: this !> zone water object
    integer(kind=4), intent(in) :: tar_wat_idx !> global index of the target water being mixed
    integer(kind=4), intent(in) :: zone_id !> reactive zone id of the sub-mixture
    real(kind=8), intent(in) :: lambda_c !> zone total concentration mixing ratio
    type(aqueous_chemistry_c), intent(in) :: aq_chem !> speciated sub-mixture
    !type(solid_chemistry_c), intent(in) :: zone_solid !> its scratch solid chemistry

    !type(zone_water_out_c), allocatable, target :: tmp(:)
    integer(kind=4) :: n

    ! !> every (target water, zone) pair is pushed once per step: seeing one again means a new step began
    ! if (num_zone_waters_out>0) then
    !     if (any(zone_waters_out(1:num_zone_waters_out)%tar_wat_idx==tar_wat_idx .and. &
    !             zone_waters_out(1:num_zone_waters_out)%zone_id==zone_id)) num_zone_waters_out=0
    ! end if
    ! if (.not. allocated(zone_waters_out)) allocate(zone_waters_out(8))
    ! if (num_zone_waters_out==size(zone_waters_out)) then
    !     allocate(tmp(2*size(zone_waters_out)))
    !     tmp(1:num_zone_waters_out)=zone_waters_out(1:num_zone_waters_out)
    !     call move_alloc(tmp,zone_waters_out)
    ! end if
    ! n=num_zone_waters_out+1
    this%tar_wat_idx=tar_wat_idx
    this%zone_id=zone_id
    this%lambda_c=lambda_c
    call this%aq_chem%copy_aq_chem(aq_chem)
    !zone_waters_out(n)%solid=zone_solid
    !> the writer reaches the zone's speciation algebra through this solid chemistry
    !if (.not. associated(zone_solid%reactive_zone)) error stop "push_zone_water_out: zone solid chemistry not associated to a reactive zone"
    !call zone_waters_out(n)%solid%set_reactive_zone(zone_solid%reactive_zone)
    !if (associated(zone_solid%mineral_zone)) call zone_waters_out(n)%solid%set_mineral_zone(zone_solid%mineral_zone)
    !num_zone_waters_out=n
end subroutine set_zone_water_attributes

!> ================================================================================================
!> CURRENT-STEP (theta / new-step) reactive-mixing contribution to r_tilde.
!> Each routine adds ONLY the theta term; the (1-theta) previous-step term is added separately by
!> add_prev_step_kin_react_mix_contrib_{Euler,Lagr}. The new step is entirely at the CURRENT reactive
!> zone, so the theta term is identical for the Eulerian and Lagrangian consistent formulations (no
!> aqueous/solid phase split) -> one routine per downstream-estimation option, shared by both.
!> Index-based (all_waters + mix_idx) to match the *_idx callers. r_tilde is updated in place.
!> ================================================================================================

!> opt1 new-step term: upstream neighbours (full new rate); each downstream water reuses the closest
!> upstream water's new rate when it has the same (full) reaction set. Mixing waters with a zero
!> reaction mixing ratio (externals of the conc-aligned stencil) are skipped.
subroutine add_curr_step_react_mix_contrib_opt1(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,r_tilde)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: mix_idx(:) !> indices of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:) !> variable activity species indices of the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:) !> reaction mixing ratios, position-aligned with mix_idx (0 for conc-only waters)
    integer(kind=4), intent(in) :: num_up !> number of upstream waters
    integer(kind=4), intent(in) :: num_down !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: r_tilde(:) !> reaction rate contributions (updated in place)
    integer(kind=4) :: j !> index mixing waters
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of the current mixing water
    real(kind=8), allocatable :: S_nc_up(:,:) !> full new stoichiometric matrix of the closest upstream water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a downstream water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of the current mixing water
    real(kind=8), allocatable :: r_new_up(:) !> full new reaction rates of the closest upstream water
    real(kind=8), allocatable :: r_old(:) !> full old reaction rates of a downstream water
    real(kind=8) :: tol !> tolerance
    logical :: same_reacts !> true iff the downstream water shares the closest upstream water's full reaction set
    logical :: have_up !> true once an upstream water with a nonzero reaction ratio has been processed
    tol=all_waters(mix_idx(1))%solid_chemistry%reactive_zone%CV_params%zero
    have_up=.false.
    do j=1,num_up
        if (mixing_ratios_reacts(j)==0d0) cycle
        S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz) !> full new stoichiometric matrix, upstream water j
        r_new=all_waters(mix_idx(j))%get_r_new() !> full new reaction rates (equilibrium + kinetic)
        r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
        have_up=.true.
    end do
    if (have_up) then
        !> closest upstream water = last upstream processed
        S_nc_up=S_nc
        r_new_up=r_new
    end if
    do j=1,num_down
        if (mixing_ratios_reacts(num_up+j)==0d0) cycle
        same_reacts=.false.
        if (have_up) then
            S_nc=all_waters(mix_idx(num_up+j))%get_S_nc(ind_var_act_species_rz) !> full stoichiometric matrix of the downstream water
            if (all(shape(S_nc)==shape(S_nc_up))) then !> same reaction set requires identical full stoichiometry
                if (norm_mat_inf(S_nc-S_nc_up)<tol) same_reacts=.true.
            end if
        end if
        if (same_reacts) then
            !> reuse the closest upstream water's new reaction rate as the downstream estimate
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_up),r_new_up)
        else
            !> no matching upstream: estimate the downstream new-step rate by its own previous-step rate
            S_nc_old=all_waters(mix_idx(num_up+j))%get_S_nc_old(ind_var_act_species_rz)
            r_old=all_waters(mix_idx(num_up+j))%get_r_old()
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_old),r_old)
        end if
    end do
end subroutine

!> opt2 new-step term: upstream neighbours (full new rate) + downstream first-derivative extrapolation
!> r_est=(1+mu) r_old - mu r_old_old (mu=est_prm), transformed by the downstream old stoichiometry.
!> Mixing waters with a zero reaction mixing ratio are skipped.
subroutine add_curr_step_react_mix_contrib_opt2(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,r_tilde)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: mix_idx(:) !> indices of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:) !> variable activity species indices of the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:) !> reaction mixing ratios, position-aligned with mix_idx (0 for conc-only waters)
    integer(kind=4), intent(in) :: num_up !> number of upstream waters
    integer(kind=4), intent(in) :: num_down !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: r_tilde(:) !> reaction rate contributions (updated in place)
    integer(kind=4) :: j !> index mixing waters
    real(kind=8), allocatable :: r_est(:) !> estimated full k+1 reaction rates of a downstream water
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a downstream water
    do j=1,num_up
        if (mixing_ratios_reacts(j)==0d0) cycle
        S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz)
        r_new=all_waters(mix_idx(j))%get_r_new()
        r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
    end do
    do j=1,num_down
        if (mixing_ratios_reacts(num_up+j)==0d0) cycle
        S_nc_old=all_waters(mix_idx(num_up+j))%get_S_nc_old(ind_var_act_species_rz)
        r_est=(1d0+all_waters(mix_idx(num_up+j))%solid_chemistry%reactive_zone%CV_params%est_prm)*&
            all_waters(mix_idx(num_up+j))%get_r_old()-&
            all_waters(mix_idx(num_up+j))%solid_chemistry%reactive_zone%CV_params%est_prm*&
            all_waters(mix_idx(num_up+j))%get_r_old_old()
        r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_old),r_est)
    end do
end subroutine

!> opt3 new-step term: upstream neighbours (full new rate) + downstream second-derivative extrapolation
!> r_est=3 r_old-3 r_old_old+r_old_old_old, transformed by the downstream old stoichiometry.
!> Mixing waters with a zero reaction mixing ratio are skipped.
subroutine add_curr_step_react_mix_contrib_opt3(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,r_tilde)
    implicit none
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array
    integer(kind=4), intent(in) :: mix_idx(:) !> indices of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:) !> variable activity species indices of the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:) !> reaction mixing ratios, position-aligned with mix_idx (0 for conc-only waters)
    integer(kind=4), intent(in) :: num_up !> number of upstream waters
    integer(kind=4), intent(in) :: num_down !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: r_tilde(:) !> reaction rate contributions (updated in place)
    integer(kind=4) :: j !> index mixing waters
    real(kind=8), allocatable :: r_est(:) !> estimated full k+1 reaction rates of a downstream water
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a downstream water
    do j=1,num_up
        if (mixing_ratios_reacts(j)==0d0) cycle
        S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz)
        r_new=all_waters(mix_idx(j))%get_r_new()
        r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
    end do
    do j=1,num_down
        if (mixing_ratios_reacts(num_up+j)==0d0) cycle
        S_nc_old=all_waters(mix_idx(num_up+j))%get_S_nc_old(ind_var_act_species_rz)
        r_est=3d0*all_waters(mix_idx(num_up+j))%get_r_old()-3d0*all_waters(mix_idx(num_up+j))%get_r_old_old()+&
            all_waters(mix_idx(num_up+j))%get_r_old_old_old()
        r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_old),r_est)
    end do
end subroutine

!> This subroutine computes the reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters, for the first time step
!> The estimation of the downstream waters reaction rates is done by taking the reaction rate of the closest upstream water
subroutine compute_r_tilde_impl_opt1(mix_waters,ind_var_act_species_rz,&
                                     mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde,curr_water,lambda_self)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> target waters that mix with each other, ordered in upstream & downstream waters
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> same order as "mix_waters"
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    type(aqueous_chemistry_c), intent(in) :: curr_water !> current spatial cell (unused in this option; required for interface compatibility)
    real(kind=8), intent(in) :: lambda_self !> self reaction mixing ratio (unused in this option; required for interface compatibility)
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: S_nc(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the current mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full (equilibrium + kinetic) old stoichiometric matrix of the current mixing water
    real(kind=8), allocatable :: S_nc_up(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the closest upstream water
    real(kind=8), allocatable :: S_nc_up_old(:,:) !> full (equilibrium + kinetic) old stoichiometric matrix of the closest upstream water
    real(kind=8), allocatable :: r_new(:) !> full (equilibrium + kinetic) new reaction rates of the current mixing water
    real(kind=8), allocatable :: r_old(:) !> full (equilibrium + kinetic) old reaction rates of the current mixing water
    real(kind=8), allocatable :: r_new_up(:) !> full (equilibrium + kinetic) new reaction rates of the closest upstream water
    real(kind=8), allocatable :: r_old_up(:) !> full (equilibrium + kinetic) old reaction rates of the closest upstream water
    real(kind=8) :: tol !> tolerance
    logical :: same_reacts !> true iff the downstream water has the same reaction set (full stoichiometry) as the closest upstream water
    
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt1"
    
!> Process
    tol=mix_waters(1)%solid_chemistry%reactive_zone%CV_params%zero
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> Upstream waters: full (equilibrium + kinetic) new & old reaction rates
        do j=1,num_up
            S_nc=mix_waters(j)%get_S_nc(ind_var_act_species_rz) !> full new stoichiometric matrix, upstream water j
            S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> full old stoichiometric matrix, upstream water j
            r_new=mix_waters(j)%get_r_new() !> full new reaction rates (equilibrium + kinetic)
            r_old=mix_waters(j)%get_r_old() !> full old reaction rates (equilibrium + kinetic)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*(theta*matmul(transpose(S_nc),r_new)+&
                (1d0-theta)*matmul(transpose(S_nc_old),r_old))
        end do
    !> Downstream waters
        if (num_up>0) then
            !> closest upstream water = last upstream processed; keep its full stoichiometry and rates
            S_nc_up=S_nc
            S_nc_up_old=S_nc_old
            r_new_up=r_new
            r_old_up=r_old
            do j=1,num_down
                S_nc=mix_waters(num_up+j)%get_S_nc(ind_var_act_species_rz) !> full stoichiometric matrix of the downstream water
                same_reacts=.false.
                if (all(shape(S_nc)==shape(S_nc_up))) then !> same reaction set requires identical full stoichiometry (equilibrium + kinetic)
                    if (norm_mat_inf(S_nc-S_nc_up)<tol) same_reacts=.true.
                end if
                if (same_reacts) then
                    !> same reaction set: reuse the closest upstream water's full reaction rate
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*(theta*matmul(transpose(S_nc_up),r_new_up)+&
                        (1d0-theta)*matmul(transpose(S_nc_up_old),r_old_up))
                else
                    !> different reaction set: use the downstream water's own previous-time-step full rate
                    r_old=mix_waters(num_up+j)%get_r_old()
                    S_nc_old=mix_waters(num_up+j)%get_S_nc_old(ind_var_act_species_rz)
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*matmul(transpose(S_nc_old),r_old)
                end if
            end do
        else
            do j=1,num_down
                r_old=mix_waters(j)%get_r_old() !> full previous-time-step reaction rates of the downstream water
                S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> full old stoichiometric matrix of the downstream water
                r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(S_nc_old),r_old)
            end do
        end if
end subroutine

!> This subroutine computes the kinetic reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by extrapolation with the two previous reaction rates
subroutine compute_r_tilde_impl_opt2(mix_waters,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde,curr_water,lambda_self)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    type(aqueous_chemistry_c), intent(in) :: curr_water !> current (self) cell; contributes its (1-theta) old reaction rate to r_tilde (NOT part of mix_waters)
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; weights its (1-theta) old reaction rate
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: r_est(:) !> estimated FULL (equilibrium + kinetic) k+1 reaction rates of a downstream water (option 2: first-derivative extrapolation)
    real(kind=8), allocatable :: S_nc(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a mixing water (previous-step reactive zone)
    real(kind=8), allocatable :: r_old(:) !> full old reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_c_old(:,:) !> full old stoichiometric matrix of the current (self) cell
    real(kind=8), allocatable :: r_old_c(:) !> full old reaction rates of the current (self) cell
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    !real(kind=8) :: mu=1d0 !> parameter for extrapolation of downstream waters rk
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt2"
    
!> Process
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    !rk_tilde_down=0d0 !> initialise kinetic reaction contributions after mixing
    r_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Step (b): NEW-STEP (theta), upstream neighbours, full new rates (by reactive zone).
        do j=1,num_up
            S_nc=mix_waters(j)%get_S_nc(ind_var_act_species_rz) !> new stoichiometric matrix, upstream water j
            r_new=mix_waters(j)%get_r_new() !> new reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
        end do
    !> Step (a): PREVIOUS-STEP (1-theta), mixing neighbours (upstream + downstream) only, full old
    !> rate. The current cell (self) is excluded: its previous-step term is added unscaled by the
    !> Newton source_old, so it stays out of this Patankar-limited r_tilde.
        do j=1,num_up+num_down
            S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> old stoichiometric matrix, mixing water j
            r_old=mix_waters(j)%get_r_old() !> old reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*(1d0-theta)*matmul(transpose(S_nc_old),r_old)
        end do
    !> Step (c): NEW-STEP (theta) DOWNSTREAM estimation (option 2: first-derivative extrapolation),
    !> applied to the FULL reaction rate (equilibrium AND kinetic) by reactive zone. Each downstream
    !> water's (1-theta) old contribution is already in step (a); here we add only the theta part of
    !> the estimated rate r_est = (1+mu) r^k - mu r^{k-1} (mu = est_prm): kinetic entries use
    !> rk_old/rk_old_old, equilibrium entries use re_mean/rxn_rate_eq_old (via get_r_old/get_r_old_old).
        do j=1,num_down
            S_nc_old=mix_waters(num_up+j)%get_S_nc_old(ind_var_act_species_rz) !> stoichiometric matrix, downstream water j
            r_est=(1d0+mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%est_prm)*&
                mix_waters(num_up+j)%get_r_old()-&
                mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%est_prm*&
                mix_waters(num_up+j)%get_r_old_old()
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_old),r_est)
        end do
end subroutine 

!> Option 4 new-step reaction-rate mixing term for the consistent WMA formulations (Eulerian and
!> Lagrangian share it: any difference is carried by the mixing ratios passed in). r_tilde is the
!> theta_r (new-step) term ONLY, built from the upstream neighbours' FULL new reaction rates; no
!> downstream new-step estimate (downstream cells enter only through the previous-step term that
!> the consistent mixers already fold into c_mix).
!>
!> Index-based: receives the full waters array (all_waters) plus an index list (mix_idx) of the
!> mixing waters instead of an already-gathered mix_waters(:) section, avoiding the expensive
!> vector-subscript gather that deep-copies ~25 allocatable components per mixing water on every
!> cell and time step.
subroutine compute_r_tilde_impl_opt4_idx(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta_r,r_tilde,curr_water_idx,lambda_self,wat_indices_zones)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array; the mixing waters are all_waters(mix_idx(:)), ordered upstream then downstream. Passing the whole array plus an index list (instead of a gathered mix_waters section) avoids a costly vector-subscript gather that deep-copies ~25 allocatable components per mixing water every cell/time step.
    integer(kind=4), intent(in) :: mix_idx(:) !> indices into all_waters of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters (opt4 has no downstream new-step term)
    real(kind=8), intent(in) :: theta_r !> time weighting factor of the KINETIC reactions (equilibrium block enters with weight 1)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    integer(kind=4), intent(in) :: curr_water_idx !> index into all_waters of the current (self) cell; UNUSED (its previous-step term is folded into c_mix by the mixer), kept for interface compatibility with the other *_idx variants
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; UNUSED (see curr_water_idx), kept for interface compatibility with the other *_idx variants
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (cols): global water indices in that zone (col_1); UNUSED here, kept for interface compatibility
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_eq !> number of equilibrium reactions (leading block of S_nc / r_new)
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of an upstream mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of an upstream mixing water
    
    num_mix_wat=size(mix_idx)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt4_idx"
    
!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> --------------------------------------------------------------------------------------------
    !> Consistent formulation (Eulerian and Lagrangian alike: the difference is carried by the mixing
    !> ratios). r_tilde is the NEW-STEP (theta_r) term ONLY, from the upstream neighbours' FULL new
    !> reaction rates (equilibrium + kinetic; no downstream new-step term). The previous-step
    !> (1-theta_r) contribution is already folded into c_mix by the consistent mixers, so adding it
    !> here would double-count it. Mixing waters with a zero reaction mixing ratio (the externals of
    !> the conc-aligned stencil) are skipped.
    !> --------------------------------------------------------------------------------------------
    do j=1,num_up
        if (mixing_ratios_reacts(j)==0d0) cycle
        S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz) !> full new stoichiometric matrix, upstream water j
        r_new=all_waters(mix_idx(j))%get_r_new() !> full new reaction rates (equilibrium + kinetic)
        n_eq=all_waters(mix_idx(j))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        r_new(n_eq+1:)=theta_r*r_new(n_eq+1:) !> theta_r weights only the kinetic tail
        r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(S_nc),r_new)
    end do
end subroutine

!> Index-based variant of compute_r_tilde_impl_opt1 (closest-upstream downstream estimate): receives
!> the full waters array plus an index list of the mixing waters instead of a gathered mix_waters(:)
!> section, avoiding the costly vector-subscript gather. curr_water/lambda_self are unused here
!> (kept for interface compatibility with the other *_idx variants).
subroutine compute_r_tilde_impl_opt1_idx(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,r_tilde,curr_water_idx,lambda_self,wat_indices_zones)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array; the mixing waters are all_waters(mix_idx(:)), ordered upstream then downstream
    integer(kind=4), intent(in) :: mix_idx(:) !> indices into all_waters of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    integer(kind=4), intent(in) :: curr_water_idx !> index into all_waters of the current (self) cell; UNUSED (its previous-step term is folded into c_mix by the mixer), kept for interface compatibility with the other *_idx variants
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; UNUSED (see curr_water_idx), kept for interface compatibility with the other *_idx variants
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (cols): global target-water indices in that zone (col_1); used for the previous-step by-zone term
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: S_nc(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the current mixing water
    real(kind=8), allocatable :: S_nc_up(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the closest upstream water
    real(kind=8), allocatable :: r_new(:) !> full (equilibrium + kinetic) new reaction rates of the current mixing water
    real(kind=8), allocatable :: r_new_up(:) !> full (equilibrium + kinetic) new reaction rates of the closest upstream water
    real(kind=8) :: tol !> tolerance
    logical :: same_reacts !> true iff the downstream water has the same reaction set (full stoichiometry) as the closest upstream water

    num_mix_wat=size(mix_idx)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt1_idx"

!> Process
    tol=all_waters(mix_idx(1))%solid_chemistry%reactive_zone%CV_params%zero
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> NEW-STEP (theta) term ONLY (upstream full new rate; downstream reuses closest-upstream new rate).
    !> The previous-step (1-theta) term of self and neighbours is already folded into c_mix by the
    !> consistent mixers (compute_c_mix_by_zones_cons / compute_c_mix_global_cons); adding it here too
    !> would double-count it in c_hat = c_mix + Delta_t*r_tilde.
    call add_curr_step_react_mix_contrib_opt1(all_waters,mix_idx,ind_var_act_species_rz,&
        mixing_ratios_reacts,num_up,num_down,theta,r_tilde)
end subroutine

!> Index-based variant of compute_r_tilde_impl_opt2 (first-derivative downstream extrapolation):
!> receives the full waters array plus an index list of the mixing waters instead of a gathered
!> mix_waters(:) section, avoiding the costly vector-subscript gather.
subroutine compute_r_tilde_impl_opt2_idx(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta_r,r_tilde,curr_water_idx,lambda_self,wat_indices_zones)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array; the mixing waters are all_waters(mix_idx(:)), ordered upstream then downstream
    integer(kind=4), intent(in) :: mix_idx(:) !> indices into all_waters of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta_r !> time weighting factor of the KINETIC reactions (equilibrium block enters with weight 1)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    integer(kind=4), intent(in) :: curr_water_idx !> index into all_waters of the current (self) cell; UNUSED (its previous-step term is folded into c_mix by the mixer), kept for interface compatibility with the other *_idx variants
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; UNUSED (see curr_water_idx), kept for interface compatibility with the other *_idx variants
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (cols): global target-water indices in that zone (col_1); UNUSED here, kept for interface compatibility
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_eq !> number of equilibrium reactions (leading block of the rate vector / stoichiometric matrix)
    real(kind=8), allocatable :: r_est(:) !> estimated FULL k+1 reaction rates of a downstream water (first-derivative extrapolation)
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a mixing water (previous-step reactive zone)

    num_mix_wat=size(mix_idx)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt2_idx"

!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> r_tilde is the NEW-STEP term ONLY: the previous-step (1-theta_r) term of self and neighbours
    !> is already folded into c_mix by the consistent mixers, so it must not be added here.
    !> Mixing waters with a zero reaction mixing ratio (externals of the conc-aligned stencil) are skipped.
    !> Step (b): NEW-STEP, upstream neighbours, full new rates (eq weight 1, kin weight theta_r).
        do j=1,num_up
            if (mixing_ratios_reacts(j)==0d0) cycle
            S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz)
            r_new=all_waters(mix_idx(j))%get_r_new()
            n_eq=all_waters(mix_idx(j))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
            r_new(n_eq+1:)=theta_r*r_new(n_eq+1:) !> theta_r weights only the kinetic tail
            r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(S_nc),r_new)
        end do
    !> Step (c): NEW-STEP DOWNSTREAM first-derivative estimation r_est=(1+mu) r_old - mu r_old_old (mu=est_prm),
    !> in the downstream water's previous-step layout (r_old/r_old_old/S_nc_old rows: get_r_hist / get_S_nc_old).
        do j=1,num_down
            if (mixing_ratios_reacts(num_up+j)==0d0) cycle
            S_nc_old=all_waters(mix_idx(num_up+j))%get_S_nc_old(ind_var_act_species_rz)
            r_est=(1d0+all_waters(mix_idx(num_up+j))%solid_chemistry%reactive_zone%CV_params%est_prm)*&
                all_waters(mix_idx(num_up+j))%get_r_old()-&
                all_waters(mix_idx(num_up+j))%solid_chemistry%reactive_zone%CV_params%est_prm*&
                all_waters(mix_idx(num_up+j))%get_r_old_old()
            n_eq=size(r_est)-all_waters(mix_idx(num_up+j))%indices_rk_old%num_cols !> equilibrium rows precede the kinetic tail
            r_est(n_eq+1:)=theta_r*r_est(n_eq+1:) !> theta_r weights only the kinetic tail
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*matmul(transpose(S_nc_old),r_est)
        end do
end subroutine

!> Index-based variant of compute_r_tilde_impl_opt3 (second-derivative downstream extrapolation):
!> receives the full waters array plus an index list of the mixing waters instead of a gathered
!> mix_waters(:) section, avoiding the costly vector-subscript gather.
subroutine compute_r_tilde_impl_opt3_idx(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta_r,r_tilde,curr_water_idx,lambda_self,wat_indices_zones)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array; the mixing waters are all_waters(mix_idx(:)), ordered upstream then downstream
    integer(kind=4), intent(in) :: mix_idx(:) !> indices into all_waters of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta_r !> time weighting factor of the KINETIC reactions (equilibrium block enters with weight 1)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    integer(kind=4), intent(in) :: curr_water_idx !> index into all_waters of the current (self) cell; UNUSED (its previous-step term is folded into c_mix by the mixer), kept for interface compatibility with the other *_idx variants
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; UNUSED (see curr_water_idx), kept for interface compatibility with the other *_idx variants
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (cols): global target-water indices in that zone (col_1); UNUSED here, kept for interface compatibility
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_eq !> number of equilibrium reactions (leading block of the rate vector / stoichiometric matrix)
    real(kind=8), allocatable :: r_est(:) !> estimated FULL k+1 reaction rates of a downstream water (second-derivative extrapolation)
    real(kind=8), allocatable :: S_nc(:,:) !> full new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a mixing water (previous-step reactive zone)

    num_mix_wat=size(mix_idx)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt3_idx"

!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> r_tilde is the NEW-STEP term ONLY: the previous-step (1-theta_r) term of self and neighbours
    !> is already folded into c_mix by the consistent mixers, so it must not be added here.
    !> Mixing waters with a zero reaction mixing ratio (externals of the conc-aligned stencil) are skipped.
    !> Step (b): NEW-STEP, upstream neighbours, full new rates (eq weight 1, kin weight theta_r).
        do j=1,num_up
            if (mixing_ratios_reacts(j)==0d0) cycle
            S_nc=all_waters(mix_idx(j))%get_S_nc(ind_var_act_species_rz)
            r_new=all_waters(mix_idx(j))%get_r_new()
            n_eq=all_waters(mix_idx(j))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
            r_new(n_eq+1:)=theta_r*r_new(n_eq+1:) !> theta_r weights only the kinetic tail
            r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(S_nc),r_new)
        end do
    !> Step (c): NEW-STEP DOWNSTREAM second-derivative estimation r_est=3 r_old-3 r_old_old+r_old_old_old,
    !> in the downstream water's previous-step layout (r_old*/S_nc_old rows: get_r_hist / get_S_nc_old).
        do j=1,num_down
            if (mixing_ratios_reacts(num_up+j)==0d0) cycle
            S_nc_old=all_waters(mix_idx(num_up+j))%get_S_nc_old(ind_var_act_species_rz)
            r_est=3d0*all_waters(mix_idx(num_up+j))%get_r_old()-3d0*all_waters(mix_idx(num_up+j))%get_r_old_old()+&
                all_waters(mix_idx(num_up+j))%get_r_old_old_old()
            n_eq=size(r_est)-all_waters(mix_idx(num_up+j))%indices_rk_old%num_cols !> equilibrium rows precede the kinetic tail
            r_est(n_eq+1:)=theta_r*r_est(n_eq+1:) !> theta_r weights only the kinetic tail
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*matmul(transpose(S_nc_old),r_est)
        end do
end subroutine

!> Lagrangian consistent variant of compute_r_tilde_impl_opt1_idx (closest-upstream downstream
!> estimate). The previous-step term uses the same reactive-zone previous-step helper as the Eulerian
!> version (no aqueous/solid phase split), so this routine is independent of the Euler vs Lagrangian
!> distinction: any difference is carried entirely by the mixing ratios passed in. The theta
!> (new-step) term is at the current reactive zone.
subroutine compute_r_tilde_impl_opt1_idx_Lagr(all_waters,mix_idx,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,r_tilde,curr_water_idx,lambda_self,wat_indices_zones)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: all_waters(:) !> full waters array; the mixing waters are all_waters(mix_idx(:)), ordered upstream then downstream
    integer(kind=4), intent(in) :: mix_idx(:) !> indices into all_waters of the mixing waters (upstream then downstream)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    integer(kind=4), intent(in) :: curr_water_idx !> index into all_waters of the current (self) cell; UNUSED (its previous-step term is folded into c_mix by the mixer), kept for interface compatibility with the other *_idx variants
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; UNUSED (see curr_water_idx), kept for interface compatibility with the other *_idx variants
    type(int_array_c), intent(in) :: wat_indices_zones !> per reactive zone (cols): global target-water indices in that zone (col_1); used for the previous-step by-zone term
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: S_nc(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the current mixing water
    real(kind=8), allocatable :: S_nc_up(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of the closest upstream water
    real(kind=8), allocatable :: r_new(:) !> full (equilibrium + kinetic) new reaction rates of the current mixing water
    real(kind=8), allocatable :: r_new_up(:) !> full (equilibrium + kinetic) new reaction rates of the closest upstream water
    real(kind=8) :: tol !> tolerance
    logical :: same_reacts !> true iff the downstream water has the same reaction set (full stoichiometry) as the closest upstream water

    num_mix_wat=size(mix_idx)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt1_idx_Lagr"

!> Process
    tol=all_waters(mix_idx(1))%solid_chemistry%reactive_zone%CV_params%zero
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> NEW-STEP (theta) term ONLY (upstream full new rate; downstream reuses closest-upstream new rate).
    !> The previous-step (1-theta) term is already folded into c_mix by the consistent mixers.
    call add_curr_step_react_mix_contrib_opt1(all_waters,mix_idx,ind_var_act_species_rz,&
        mixing_ratios_reacts,num_up,num_down,theta,r_tilde)
end subroutine

!> This subroutine computes the reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by assuming there is a fraction of the reaction that occurs at the original target
!>
!> Gathered-array version: receives the mixing waters already gathered into mix_waters(:).
!> Retained for the callers that still use the procedure-pointer interface (Euler consistent
!> solver and solve_reactive_mixing_ideal_cons). The Lagrangian consistent solvers use the
!> faster index-based variant compute_r_tilde_impl_opt4_idx instead.
subroutine compute_r_tilde_impl_opt4(mix_waters,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde,curr_water,lambda_self)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> target waters that mix with each other (must be ordered in upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    type(aqueous_chemistry_c), intent(in) :: curr_water !> current (self) cell; its FULL old reaction rate is added to the (1-theta) by-reactive-zone term (Eulerian, in place); the theta term is unaffected
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; weights its FULL old reaction rate in the previous-time-step by-reactive-zone term
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_v !> number of variable activity species in downstream water
    integer(kind=4) :: n_aq_kin !> number of aqueous kinetic reactions in downstream water
    integer(kind=4) :: n_head !> number of leading non-aqueous (mineral/gas cst-act) reactions in a rate vector
    integer(kind=4) :: n_aq !> number of aqueous (eq + kin) reactions in a rate vector
    real(kind=8), allocatable :: S_nc_c(:,:) !> current-cell new stoichiometric matrix (for solid/gas reactions)
    real(kind=8), allocatable :: S_nc_c_old(:,:) !> current-cell old stoichiometric matrix (for solid/gas reactions)
    real(kind=8), allocatable :: r_new_c(:) !> current-cell new reaction rates (solid/gas only after masking)
    real(kind=8), allocatable :: r_old_c(:) !> current-cell old reaction rates (solid/gas only after masking)
    integer(kind=4) :: n_aq_kin_curr !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin_curr !> number of aqueous kinetic reactions in current water
    real(kind=8), allocatable :: r_old(:) !> old kinetic reaction rates
    real(kind=8), allocatable :: r_new(:) !> new kinetic reaction rates
    real(kind=8), allocatable :: S_nc_old(:,:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: S_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_curr(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: A(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: b(:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: conc_nc(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: conc_nc_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: conc_nc_est(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8) :: sum_lambda_up !> sum of upstream reaction mixing ratios (weight of the factored current-cell solid/gas contribution)
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt4"
    
!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> --------------------------------------------------------------------------------------------
    !> Eulerian consistent formulation.
    !> NEW-STEP (theta) term: KEPT as the neighbour (upstream) sequence with FULL new reaction rates
    !>   (no aqueous/solid masking).
    !> PREVIOUS-STEP (1-theta) term: computed BY REACTIVE ZONES over the current cell (self, when
    !>   supplied via curr_water + lambda_self) and all mixing neighbours (upstream + downstream),
    !>   each contributing its FULL old reaction rate (aqueous AND solid/gas) transformed by its own
    !>   old reactive-zone stoichiometry and weighted by its reaction mixing ratio. Because waters
    !>   sharing an old reactive zone have the same S_nc_old, this per-water sum equals the
    !>   zone-grouped form  sum_z S_z,old^T ( sum_{j in z} lambda_j r_old_j ).
    !> --------------------------------------------------------------------------------------------
        !> NEW-STEP (theta): upstream neighbours, full new rates.
        do j=1,num_up
            S_nc=mix_waters(j)%get_S_nc(ind_var_act_species_rz) !> new stoichiometric matrix, upstream water j
            r_new=mix_waters(j)%get_r_new() !> new reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
        end do
        !> PREVIOUS-STEP (1-theta): mixing neighbours (upstream + downstream) only, full old rate. The
        !> current cell (self) is excluded: its previous-step term is added unscaled by the Newton
        !> source_old, so it stays out of this Patankar-limited r_tilde.
        do j=1,num_up+num_down
            S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> old stoichiometric matrix, mixing water j
            r_old=mix_waters(j)%get_r_old() !> old reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*(1d0-theta)*matmul(transpose(S_nc_old),r_old)
        end do
end subroutine

!> This subroutine computes the kinetic reaction rate contributions after water mixing implicitly for both upstream and downstream waters
!> The estimation of the downstream waters reaction amoutns is done by extrapolation with the second time derivative of the reaction rates
!> We assume uniform time step
subroutine compute_r_tilde_impl_opt3(mix_waters,ind_var_act_species_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde,curr_water,lambda_self)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other (must be ordered in upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    type(aqueous_chemistry_c), intent(in) :: curr_water !> current (self) cell; contributes its (1-theta) old reaction rate to r_tilde (NOT part of mix_waters)
    real(kind=8), intent(in) :: lambda_self !> reaction mixing ratio of the current (self) cell; weights its (1-theta) old reaction rate
!> Variables
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: r_est(:) !> estimated FULL (equilibrium + kinetic) k+1 reaction rates of a downstream water (option 3: second-derivative extrapolation)
    real(kind=8), allocatable :: S_nc(:,:) !> full (equilibrium + kinetic) new stoichiometric matrix of a mixing water
    real(kind=8), allocatable :: r_new(:) !> full new reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_old(:,:) !> full old stoichiometric matrix of a mixing water (previous-step reactive zone)
    real(kind=8), allocatable :: r_old(:) !> full old reaction rates of a mixing water
    real(kind=8), allocatable :: S_nc_c_old(:,:) !> full old stoichiometric matrix of the current (self) cell
    real(kind=8), allocatable :: r_old_c(:) !> full old reaction rates of the current (self) cell
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt3"
    
!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> Step (b): NEW-STEP (theta), upstream neighbours, full new rates (by reactive zone).
        do j=1,num_up
            S_nc=mix_waters(j)%get_S_nc(ind_var_act_species_rz) !> new stoichiometric matrix, upstream water j
            r_new=mix_waters(j)%get_r_new() !> new reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*theta*matmul(transpose(S_nc),r_new)
        end do
    !> Step (a): PREVIOUS-STEP (1-theta), mixing neighbours (upstream + downstream) only, full old
    !> rate. The current cell (self) is excluded: its previous-step term is added unscaled by the
    !> Newton source_old, so it stays out of this Patankar-limited r_tilde.
        do j=1,num_up+num_down
            S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> old stoichiometric matrix, mixing water j
            r_old=mix_waters(j)%get_r_old() !> old reaction rates (FULL)
            r_tilde=r_tilde+mixing_ratios_reacts(j)*(1d0-theta)*matmul(transpose(S_nc_old),r_old)
        end do
    !> Step (c): NEW-STEP (theta) DOWNSTREAM estimation (option 3: second-derivative extrapolation),
    !> applied to the FULL reaction rate (equilibrium AND kinetic) by reactive zone. Each downstream
    !> water's (1-theta) old contribution is already in step (a); here we add only the theta part of the
    !> estimated rate. The quadratic (second-order) extrapolation to k+1 is r^{k+1,est}=3 r^k-3 r^{k-1}+r^{k-2};
    !> splitting off the (1-theta) old term already added leaves r_est = 3 r_old - 3 r_old_old + r_old_old_old.
        do j=1,num_down
            S_nc_old=mix_waters(num_up+j)%get_S_nc_old(ind_var_act_species_rz) !> stoichiometric matrix, downstream water j
            r_est=3d0*mix_waters(num_up+j)%get_r_old()-3d0*mix_waters(num_up+j)%get_r_old_old()+&
                mix_waters(num_up+j)%get_r_old_old_old()
            r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*theta*matmul(transpose(S_nc_old),r_est)
        end do
end subroutine 

!> This subroutine computes the kinetic reaction rate contributions after water mixing, explicitly for all target waters
subroutine compute_r_tilde_expl(mix_waters,ind_var_act_species_rz,mixing_ratios_reacts,&
    num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    !integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    !integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each mixing water
!> Pre-process
    num_mix_wat=size(mix_waters)
!> Process
    r_tilde=0d0 !> initialise kinetic reaction rates contributions after mixing
    !    !> Autentica chapuza
    !        
    do j=1,num_mix_wat
        if (mix_waters(j)%indices_rk%num_cols>0) then
            Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
            if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Sk_nc),mix_waters(j)%rxn_rate_kin_old)
            end if
            if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Sk_nc),mix_waters(j)%solid_chemistry%rxn_rate_kin_old)
            end if
            deallocate(Sk_nc)
        end if
    end do

    !!do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !!> Mixing waters
    !    !!> Upstream waters
        !    !> Downstream waters
    !!r_tilde_aux=r_tilde !> we store the aqueous kinetic reaction rates after mixing (chapuza)
    !    !> we permute the indices in "perm"
end subroutine 

!> This subroutine initialises concentrations for any iterative method
!> It takes a linear combination of the concentrations in the two previous time steps
subroutine initialise_iterative_method(conc_old_old,conc_old,param,initial_guess)
    real(kind=8), intent(in) :: conc_old_old(:) !> concentrations at time step k-1
    real(kind=8), intent(in) :: conc_old(:) !> concentrations at time step k
    real(kind=8), intent(in) :: param !> parameter for linear combination
    real(kind=8), intent(out) :: initial_guess(:) !> initial guess concentrations at time step k+1 (must be allocated before calling)
    
    integer(kind=4) :: i !> loop index
    
    if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
    do i=1,size(conc_old)
        if ((1d0+param)*conc_old(i)-param*conc_old_old(i)<0d0) then
            initial_guess(i)=conc_old(i)
        else
            initial_guess(i)=(1d0+param)*conc_old(i)-param*conc_old_old(i)
        end if
    end do
end subroutine

!> This subroutine initialises aqueous primary concentrations for any iterative method
!> It takes a linear combination of the aqueous primary concentrations in the two previous time steps
subroutine initialise_c1_aq_iterative_method(this,param)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: param !> parameter for linear combination
    
    integer(kind=4) :: i !> loop index
    
    if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        if ((1d0+param)*this%conc_old(this%ind_var_act_species(i))-param*this%conc_old_old(this%ind_var_act_species(i))<0d0) then
            this%concentrations(this%ind_var_act_species(i))=this%conc_old(this%ind_var_act_species(i))
        else
            this%concentrations(this%ind_var_act_species(i))=(1d0+param)*this%conc_old(this%ind_var_act_species(i))-&
                param*this%conc_old_old(this%ind_var_act_species(i))
        end if
    end do
end subroutine

!> \brief Sets indices of primary, secondary, and variable activity aqueous species
!!
!! \details This subroutine establishes the indices of aqueous species in the concentrations attribute of the aqueous chemistry object.
!! It creates three key index arrays:
!! - **ind_prim_species**: Indices of primary (basis) aqueous species in concentrations array
!! - **ind_sec_species**: Indices of secondary aqueous species in concentrations array (complexes + constant activity aqueous species)
!! - **ind_var_act_species**: Indices of all aqueous variable activity species in concentrations array (primary + secondary variable)
!!
!! **Two modes of operation based on flag_comp**:
!!
!! **Case 1: flag_comp = .false. (Standard speciation algebra)**
!! - Assumes aqueous species are already sorted: primary species first, then secondary species
!! - Primary species indices: 1 to num_aq_prim_species
!! - Secondary species indices: num_aq_prim_species+1 to num_aq_species
!! - Variable activity species are identified by checking cst_act_flag for each species
!!
!! **Case 2: flag_comp = .true. (Component-based speciation algebra)**
!! - Handles constant activity species that may appear among the first num_aq_prim_species species
!! - Uses cntr_cst_act_sp to track how many constant activity species are in the "primary" positions
!! - Adjusts index assignment to ensure true primary species are correctly identified
!! - Secondary species array includes: variable activity secondary species + all constant activity species
!!
!! **Index array structure**:
!! - ind_prim_species[1:num_aq_prim_species]: Maps to actual primary species
!! - ind_sec_species[1:num_sec_aq_species]: Maps to secondary species (including constant activity)
!! - ind_var_act_species[1:num_aq_var_act_species]: Maps to all variable activity species
!!   * First num_aq_prim_species entries: primary variable activity species
!!   * Remaining entries: secondary variable activity species
!!
!! \param[in,out] this Aqueous chemistry object with species indices to be set
!!
!! \note The subroutine assumes that indices_aq_species is already set and maps aqueous species
!!       to their positions in the global chemical system species array.
!!
!! \warning Deallocates and reallocates index arrays if flag_comp = .true., so any previous
!!          index information is lost.
!!
!! \see speciation_alg, indices_aq_species
!!
    subroutine set_ind_species(this)!,ind_prim_species,ind_sec_var_act_species)
    class(aqueous_chemistry_c) :: this          !< [in,out] Aqueous chemistry object
    !integer(kind=4), intent(in), optional :: ind_prim_species(:), ind_sec_var_act_species(:) !> if one of them is present, the other must be present as well
    
!> \name Local variables
!! @{
    integer(kind=4) :: i                        !< Loop counter for all aqueous species
    integer(kind=4) :: j                        !< Counter for variable activity primary species (flag_comp=false) or all primary species (flag_comp=true)
    integer(kind=4) :: k                        !< Counter for variable activity secondary species
    integer(kind=4) :: l                        !< Counter for constant activity species (secondary species array)
    integer(kind=4) :: cntr_cst_act_sp          !< Counter for constant activity species in first num_aq_prim_species positions (flag_comp=true only)
    integer(kind=4), allocatable :: aux_ind_sp(:) !< Auxiliary index array (not currently used, kept for potential future use)
    integer :: n_p_aq_old                       !< Number of primary aqueous species before changing definition (not currently used)
!! @}
    if (allocated(this%ind_prim_species)) then
        deallocate(this%ind_prim_species)
        !deallocate(this%ind_prim_species_old)
    end if
    if (allocated(this%ind_sec_species)) then
        deallocate(this%ind_sec_species)
    end if
    if (allocated(this%ind_var_act_species)) then
        deallocate(this%ind_var_act_species)
    end if
    if (.not. associated(this%solid_chemistry)) then
        error stop "set_ind_species: solid_chemistry pointer is not associated"
    end if
    if (.not. associated(this%solid_chemistry%reactive_zone)) then
        error stop "set_ind_species: reactive_zone pointer is not associated"
    end if
    flush(6) 
    allocate(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    allocate(this%ind_prim_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    allocate(this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
        
!> \subsection init_counters Initialize counters
!! All counters are initialized to 1 (Fortran arrays are 1-indexed).
!! - j: Tracks position in ind_prim_species or ind_var_act_species
!! - k: Tracks position in secondary variable activity species
!! - l: Tracks position in constant activity species (placed at end of ind_sec_species)
        j=1                                     !< Initialize counter for primary species (or primary variable activity)
        k=1                                     !< Initialize counter for secondary variable activity species
        l=1                                     !< Initialize counter for secondary species (or constant activity species)
        
!> \subsection case_flag_comp_false Case 1: Standard speciation algebra (flag_comp = .false.)
!! In this case, aqueous species are assumed to be already sorted:
!! - First num_aq_prim_species species are primary (basis) species
!! - Remaining species are secondary (complexes)
!! Variable activity status is determined by checking cst_act_flag for each species.
        if (this%solid_chemistry%reactive_zone%speciation_alg%flag_comp .eqv. .false.) then
            !> Allocate index arrays based on speciation algebra dimensions
            
            !> \name Loop 1: Process primary species
            !! Iterate through first num_prim_species species. All are assigned to ind_prim_species.
            !! Variable activity primary species are also added to ind_var_act_species.
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                this%ind_prim_species(i)=this%indices_aq_species(i)  !< Store index in chemical system
                if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .false.) then  !< Check if variable activity
                    this%ind_var_act_species(j)=this%indices_aq_species(i)        !< Store in variable activity array
                    j=j+1                       !< Increment variable activity counter
                end if
            end do
            !> \name Loop 2: Process secondary species
            !! Iterate through remaining species (num_aq_prim_species+1 to num_aq_species).
            !! All are assigned to ind_sec_species. Variable activity secondary species are
            !! also added to ind_var_act_species (after primary variable activity species).
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                !> Store index of secondary species i (offset by num_aq_prim_species)
                this%ind_sec_species(i)=this%indices_aq_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                !> Check if this secondary species has variable activity
                if (this%aq_phase%aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%cst_act_flag &
                    .eqv. .false.) then
                    !> Store in variable activity array (continuing from where primary species left off)
                    this%ind_var_act_species(j)=this%indices_aq_species(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                    j=j+1                       !< Increment variable activity counter
                end if
            end do
            
!> \subsection case_flag_comp_true Case 2: Component-based speciation algebra (flag_comp = .true.)
!! This case handles situations where constant activity species may appear among the first
!! num_aq_prim_species positions in the species array. This requires special handling to ensure
!! that only true variable activity species are assigned to ind_prim_species.
!!
!! **Algorithm logic**:
!! - cntr_cst_act_sp tracks how many constant activity species appear in the first num_aq_prim_species positions
!! - The cutoff for primary species is dynamically adjusted: num_aq_prim_species + cntr_cst_act_sp
!! - ind_sec_species is structured as: [variable activity secondary species | constant activity species]
!! - ind_var_act_species is structured as: [primary variable activity | secondary variable activity]
!!
!! **Example**: If we have 5 species designated as "primary" but 1 is constant activity (e.g., H2O):
!! - Species 1-5 are in "primary" positions, but species 3 is H2O (constant activity)
!! - cntr_cst_act_sp = 1 after processing species 3
!! - True primary species are those with i <= 5+1 = 6 that are variable activity
!! - H2O goes to end of ind_sec_species
        else
            cntr_cst_act_sp=0                   !< Initialize counter for constant activity species in "primary" positions
            
            !> \name Loop: Process all species with activity-based classification
            !! Iterate through all aqueous species and classify them based on three criteria:
            !! 1. Constant activity goes to end of ind_sec_species
            !! 2. Variable activity in "primary" positions goes to ind_prim_species and ind_var_act_species
            !! 3. Variable activity in "secondary" positions goes to ind_sec_species and ind_var_act_species
            do i=1,this%aq_phase%num_species
                !> **Case A: Constant activity species**
                !! Place at end of ind_sec_species array (after all variable activity secondary species).
                !! If this species is in a "primary" position, increment cntr_cst_act_sp to adjust cutoff.
                if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .true.) then
                    !> Store at position: num_aq_sec_var_act_species + l (fills from end)
                    this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+l)=&
                        this%indices_aq_species(i)
                    l=l+1                       !< Increment constant activity counter
                    !> Check if this constant activity species is in a "primary" position
                    if (i<=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                        cntr_cst_act_sp=cntr_cst_act_sp+1  !< Increment constant activity in primary positions counter
                    end if
                
                !> **Case B: Variable activity species in adjusted "primary" range**
                !! Species position i is within num_aq_prim_species + cntr_cst_act_sp.
                !! These are true primary species (variable activity) and go to ind_prim_species.
                else if (i<=&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+cntr_cst_act_sp) then
                    this%ind_prim_species(j)=this%indices_aq_species(i)        !< Store in primary species array
                    this%ind_var_act_species(j)=this%indices_aq_species(i)     !< Also store in variable activity array
                    j=j+1                       !< Increment primary/variable activity counter
                
                !> **Case C: Variable activity secondary species**
                !! Species beyond the adjusted primary cutoff are secondary variable activity species.
                else
                    this%ind_sec_species(k)=this%indices_aq_species(i)         !< Store in secondary species array
                    !> Also store in variable activity array (after primary species)
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+k)=&
                        this%indices_aq_species(i)
                    k=k+1                       !< Increment secondary variable activity counter
                end if
            end do
        end if  !< End of flag_comp conditional
        
!> \subsection postconditions Post-conditions after execution
!! After this subroutine completes, the following arrays are populated:
!! - **ind_prim_species**: Maps primary species numbers (1:num_aq_prim_species) to chemical system indices
!! - **ind_sec_species**: Maps secondary species numbers (1:num_sec_aq_species) to chemical system indices
!! - **ind_var_act_species**: Maps variable activity species (1:num_aq_var_act_species) to chemical system indices
!!
!! **Index array dimensions**:
!! - size(ind_prim_species) = num_aq_prim_species
!! - size(ind_sec_species) = num_sec_aq_species
!! - size(ind_var_act_species) = num_aq_var_act_species = num_aq_prim_species + num_aq_sec_var_act_species
!!
!! These index arrays are critical for:
!! - Mapping between local aqueous chemistry arrays and global chemical system arrays
!! - Identifying which species participate in mass balance equations
!! - Separating known (primary) from unknown (secondary) species in speciation calculations
    end subroutine

!> \brief set_ind_species followed by the reactive zone's speciation column permutation.
!! \details compute_speciation_alg_arrays may swap a primary and a secondary variable-activity
!! species (e.g. so4-2 <-> co3-2 in an anhydrite zone) so the eliminated-secondary block Se_nc_2
!! is non-singular; the permutation is recorded in reactive_zone%ind_var_act_species (identity
!! when no swap). set_ind_species builds the index arrays in natural aqueous order, so a water
!! that only calls it has a primary/secondary partition inconsistent with its zone's comp_mat and
!! inv_Se_nc_2, and the mass-action law gets applied to the wrong species. This routine reorders
!! the freshly built arrays with that permutation (same result as the swap read_tar_wat_line
!! applies at read time). Use it whenever the index arrays are rebuilt on an already set-up zone.
    subroutine set_ind_species_react_zone(this)
    class(aqueous_chemistry_c) :: this          !< [in,out] Aqueous chemistry object (solid_chemistry%reactive_zone associated)
    integer(kind=4) :: j                        !< speciation column position
    integer(kind=4) :: n_aq_v                   !< number of aqueous variable activity species
    integer(kind=4) :: n_p_aq                   !< number of aqueous primary species
    integer(kind=4) :: n2v_aq                   !< number of aqueous secondary variable activity species
    integer(kind=4), allocatable :: ind_nat(:)  !< ind_var_act_species in natural aqueous order (set_ind_species output)
    call this%set_ind_species()
    associate(rz => this%solid_chemistry%reactive_zone)
    if (allocated(rz%ind_var_act_species)) then
        n_aq_v=rz%speciation_alg%num_aq_var_act_species
        n_p_aq=rz%speciation_alg%num_aq_prim_species
        n2v_aq=rz%speciation_alg%num_aq_sec_var_act_species
        ind_nat=this%ind_var_act_species
        !> zone column j holds natural position rz%ind_var_act_species(j) (only aqueous columns apply here)
        do j=1,min(n_aq_v,size(rz%ind_var_act_species))
            if (rz%ind_var_act_species(j)>=1 .and. rz%ind_var_act_species(j)<=n_aq_v) then
                this%ind_var_act_species(j)=ind_nat(rz%ind_var_act_species(j))
            end if
        end do
        this%ind_prim_species=this%ind_var_act_species(1:n_p_aq)
        this%ind_sec_species(1:n2v_aq)=this%ind_var_act_species(n_p_aq+1:n_p_aq+n2v_aq)
    end if
    end associate
    end subroutine

    !> \brief Sets concentration of secondary aqueous variable activity species
    !!
    !! @param c2_aq_nc concentrations of secondary aqueous variable activity species
    !!
    subroutine set_conc_sec_aq_var_act_species(this,c2_aq_nc) !> 
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: c2_aq_nc(:)
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            this%concentrations(this%ind_var_act_species(&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))=&
            c2_aq_nc(i)
        end do
    end subroutine

    !> \brief Sets indices in the chemical system for kinetic reaction species
    !!
    !! \details Populates the indices_rk array mapping kinetic reaction species
    !! to their positions in the chemical system arrays.
    !!
    subroutine set_indices_rk(this) !> sets indices in chemical system for kinetic reactions
        !! OJO: FALTAN REACCIONES CINÃƒâ€°TICAS LINEALES
        class(aqueous_chemistry_c) :: this
        integer(kind=4) :: i,j,k,num_rk
        num_rk=0
        call this%indices_rk%allocate_array(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+&
            this%solid_chemistry%mineral_zone%num_minerals_kin+this%solid_chemistry%reactive_zone%chem_syst%num_gas_kin_reacts)
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        if (this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts>0) then !> chapuza
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
                call this%indices_rk%cols(num_rk+i)%allocate_vector(&
                    this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%num_aq_rk)
                j=1 !> counter for aqueous variable activity species
                k=1 !> counter for relevant aqueous species in kinetic reaction rates
                do
                    if (this%solid_chemistry%reactive_zone%chem_syst%species(j)%name == &
                        this%aq_phase%aq_species(&
                        this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase(k))%name) then
                        this%indices_rk%cols(num_rk+i)%col_1(k)=j
                        if (k<this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%num_aq_rk) then
                            k=k+1
                            j=1
                        else
                            exit
                        end if
                    else if (j<this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_species) then
                        j=j+1
                    else
                        error stop "Error in set_indices_rk"
                    end if
                end do
            end do
        end if
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        if (this%solid_chemistry%mineral_zone%num_minerals_kin>0) then !> chapuza
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
                call this%indices_rk%cols(num_rk+i)%allocate_vector(&
                    this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk)
                j=1 !> counter for aqueous variable activity species
                k=1 !> counter for relevant aqueous species in kinetic reaction rates
                do
                    if (this%solid_chemistry%reactive_zone%chem_syst%species(j)%name == &
                        this%aq_phase%aq_species(&
                        this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase(k))%name) then
                        this%indices_rk%cols(num_rk+i)%col_1(k)=j
                        if (k<this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk) then
                            k=k+1
                            j=1
                        else
                            exit
                        end if
                    else if (j<this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_species) then
                        j=j+1
                    else
                        error stop "Error in set_indices_rk"
                    end if
                end do
            end do
        end if
    end subroutine

    function get_rk_new(this) result(rk_new)
    !!> This function returns the new kinetic reaction rates ordered in aqueous & mineral reactions (ESTA MAL)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: rk_new(:)
    integer(kind=4) :: i
    allocate(rk_new(this%indices_rk%num_cols))
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        rk_new(i)=this%rxn_rate_kin(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rxn_rate_kin(i)
    end do
    !> Faltan gases
    end function

    function get_r_new(this) result(r_new)
    !!> This function returns the new kinetic reaction rates ordered in aqueous & mineral reactions (ESTA MAL)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_new(:)
    integer(kind=4) :: i, n_r
    allocate(r_new(this%indices_rk%num_cols+this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    n_r=0
    do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
        r_new(n_r+i)=this%solid_chemistry%rxn_rate_eq(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%num_minerals_cst_act
    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        r_new(n_r+i)=this%gas_chemistry%rxn_rate_eq(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
        r_new(n_r+i)=this%rxn_rate_eq(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
    do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
        r_new(n_r+i)=this%solid_chemistry%rxn_rate_eq(this%solid_chemistry%reactive_zone%num_minerals_cst_act+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%num_minerals_var_act
    do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
        r_new(n_r+i)=this%solid_chemistry%rxn_rate_eq(this%solid_chemistry%reactive_zone%num_minerals+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
        r_new(n_r+i)=this%gas_chemistry%rxn_rate_eq(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        r_new(n_r+i)=this%rxn_rate_kin(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        r_new(n_r+i)=this%solid_chemistry%rxn_rate_kin(i)
    end do
    !> Faltan gases
    end function

    !> Per-species contributions to the reactive-mixing mass balance of THIS target water,
    !> volume-weighted (V=this%volume) for domain integration (Appendix A form WITH reactions):
    !>   lhs    = V*(c^{k+1}-c^k)/Delta_t              storage (LHS of the RT equation)
    !>   transp = V*(c_mix-c^k)/Delta_t                transport (conservative mixing) contribution
    !>   react  = V*(c_hat-c_mix)/Delta_t + V*lambda_r*S_nc^T*r
    !>            reaction = consistent reaction-rate mixing (c_hat-c_mix = Delta_t*r_tilde) + the
    !>            equilibrium/kinetic amount (lambda_r*S_nc^T*r). The rate r uses re_mean AS-IS for the
    !>            equilibrium reactions (compute_Re already gives the step-mean extent with the CURRENT
    !>            stoichiometry: c^{k+1}-c_hat = lambda_r*S_nc^T*Re, re = Re/Delta_t) and the
    !>            theta_r weighting theta_r*rk_new only for the kinetic tail. So for the
    !>            equilibrium block transp+react = lhs exactly (residual ~ 0). For lumped solvers
    !>            c_hat = c_mix, so the mixing term vanishes.
    !> The caller sums these over all cells and reports lhs, rhs=transp+react and residual=lhs-rhs.
    subroutine get_react_mix_mass_bal_terms(this,c_mix,c_hat,lambda_r,Delta_t,theta,lhs,transp,react)
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: c_mix(:) !> transport-only mixed variable activity concentrations
        real(kind=8), intent(in) :: c_hat(:) !> concentrations after transport AND reaction-rate mixing (pre-equilibrium)
        real(kind=8), intent(in) :: lambda_r !> self reaction mixing ratio
        real(kind=8), intent(in) :: Delta_t !> time step
        real(kind=8), intent(in) :: theta !> reaction time weighting factor (theta-method)
        real(kind=8), allocatable, intent(out) :: lhs(:),transp(:),react(:)
        real(kind=8), allocatable :: c_new(:),c_old(:),S_nc(:,:),r_new(:),sr(:)
        real(kind=8) :: V
        integer(kind=4) :: n_v,n_eq_r
        n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        c_new=this%get_conc_var_act()
        c_old=this%get_conc_old_var_act_new_rz()
        S_nc=this%get_S_nc(this%solid_chemistry%reactive_zone%ind_var_act_species)
        r_new=this%get_r_new()
        n_eq_r=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        !> theta_r-weight ONLY the kinetic tail (reactions after the n_eq_r equilibrium reactions):
        !> equilibrium re is the step-mean extent from compute_Re and is used as-is, while the
        !> kinetic reaction term is multiplied by theta_r (its previous-step (1-theta_r) part is
        !> already carried in c_mix).
        r_new(n_eq_r+1:)=theta*r_new(n_eq_r+1:)
        sr=matmul(transpose(S_nc),r_new) !> S_nc^T*r: eq uses re_mean, kin uses theta_r*rk
        V=this%volume
        allocate(lhs(n_v),transp(n_v),react(n_v))
        lhs=V*(c_new(1:n_v)-c_old(1:n_v))/Delta_t
        transp=V*(c_mix(1:n_v)-c_old(1:n_v))/Delta_t
        react=V*(c_hat(1:n_v)-c_mix(1:n_v))/Delta_t + V*lambda_r*sr(1:n_v)
        !> [MBTERMS-DBG] one-shot diagnostic for the exchange mass-balance investigation
        block
            logical, save :: mbterms_dbg_done=.false.
            integer(kind=4) :: idbg
            if ((.not. mbterms_dbg_done) .and. &
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                mbterms_dbg_done=.true.
                print '(A)', "[MBTERMS-DBG] === get_react_mix_mass_bal_terms (first exchange cell) ==="
                print '(A,I0,A,I0,A,ES12.4)', "[MBTERMS-DBG] n_v=", n_v, " n_eq_r=", n_eq_r, " V=", V
                print '(A,*(I4))', "[MBTERMS-DBG] rz%ind_var_act_species=", &
                    this%solid_chemistry%reactive_zone%ind_var_act_species
                print '(A,*(ES12.4))', "[MBTERMS-DBG] r_new(1:n_eq_r)=", r_new(1:n_eq_r)
                print '(A,*(ES12.4))', "[MBTERMS-DBG] rxn_rate_eq=", this%solid_chemistry%rxn_rate_eq
                do idbg=1,n_eq_r
                    print '(A,I0,A,*(ES12.4))', "[MBTERMS-DBG] S_nc row ", idbg, " =", S_nc(idbg,1:n_v)
                end do
                print '(A,*(ES12.4))', "[MBTERMS-DBG] sr(1:n_v)=", sr(1:n_v)
                print '(A,*(ES12.4))', "[MBTERMS-DBG] c_new-c_mix=", c_new(1:n_v)-c_mix(1:n_v)
                print '(A,*(ES12.4))', "[MBTERMS-DBG] react=", react(1:n_v)
            end if
        end block
    end subroutine

    !> \brief Brings this water into equilibrium with the equilibrium minerals/gases of its reactive zone
    !!        (PHREEQC EQUILIBRIUM_PHASES style) and returns the mole transfer per primary species.
    !!
    !! \details The transport-invariant aqueous components u = U*c_nc are held fixed while the full
    !! equilibrium of the reactive zone (all constant- and variable-activity minerals and gases) is
    !! re-imposed, letting those phases dissolve/precipitate. The quantity of each aqueous primary
    !! species that must be supplied (+) or removed (-) to reach equilibrium is the change in the
    !! aqueous elemental totals
    !! \f[
    !!   \Delta n_j = T_j^{after} - T_j^{before}, \qquad T_j = \sum_i U_{aq,ji}\, c_{aq,i},
    !! \f]
    !! where \f$U_{aq}\f$ is the constant-activity component (formula) matrix restricted to aqueous
    !! species. Aqueous complexation conserves these totals, so \f$\Delta n\f$ equals the net mole
    !! transfer coming from the mineral/gas phases (the PHREEQC per-element mole transfer).
    !!
    !! On exit the water is at equilibrium: aqueous concentrations and secondary species (and activities,
    !! in the non-ideal case) are updated by the Newton speciation; the per-phase reaction amounts are
    !! stored on solid_chemistry%rxn_amnt_eq / gas_chemistry%rxn_amnt_eq (via compute_Re); and the equilibrium mineral and
    !! gas reservoir amounts are updated by the transferred moles.
    !!
    !! @param[out] Delta_n  moles of each aqueous primary species supplied(+)/removed(-) to reach equilibrium (size num_aq_prim_species)
    !! @param[out] niter    number of Newton iterations performed by the equilibrium speciation
    !! @param[out] CV_flag  .true. if the equilibrium speciation converged (Delta_n is set to zero otherwise)
    !! @param[in]  flag_ideal .true. => ideal activities, .false. => non-ideal (activity coefficients)
    subroutine equilibrate_with_phases(this,Delta_n,niter,CV_flag,flag_ideal)
        implicit none
        class(aqueous_chemistry_c) :: this !> aqueous chemistry object (the water to equilibrate)
        real(kind=8), allocatable, intent(out) :: Delta_n(:) !> per aqueous primary species mole transfer needed
        integer(kind=4), intent(out) :: niter !> Newton iterations performed
        logical, intent(out) :: CV_flag !> convergence flag
        logical, intent(in) :: flag_ideal !> ideal (.true.) vs non-ideal (.false.)
    !> Local variables
        integer(kind=4) :: i !> loop index
        integer(kind=4) :: n_p !> number of primary species (aqueous + solid + gas)
        integer(kind=4) :: n_p_aq !> number of aqueous primary species
        integer(kind=4) :: n_v !> number of variable activity species
        integer(kind=4) :: n_aq_sp !> number of aqueous species (primary + secondary)
        integer(kind=4) :: n_sec_aq !> number of secondary aqueous species
        integer(kind=4) :: n_min_kin !> number of kinetic minerals (reservoir offset)
        real(kind=8), allocatable :: U_aq(:,:) !> aqueous elemental-total (formula) matrix
        real(kind=8), allocatable :: c_aq(:) !> aqueous species concentrations (primary then secondary)
        real(kind=8), allocatable :: T_before(:),T_after(:) !> aqueous elemental totals before/after equilibration
        real(kind=8), allocatable :: c_before(:) !> variable activity species concentrations before equilibration
        real(kind=8), allocatable :: conc_nc(:) !> variable activity species concentrations at equilibrium (Newton output)
        real(kind=8), allocatable :: conc_comp(:) !> conserved (transport-invariant) components, held fixed
        real(kind=8), allocatable :: c1_ig(:) !> initial guess for primary species
        real(kind=8), allocatable :: c2v_ig(:) !> initial guess for secondary variable activity species (non-ideal path)

        n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        n_aq_sp=this%aq_phase%num_species
        n_sec_aq=n_aq_sp-n_p_aq

        !> Aqueous elemental-total operator U_aq (aq primary species x aq species), taken from the
        !> constant-activity component matrix (identity on primaries, transpose(Se_1_star) on the
        !> aqueous secondary block) exactly as the CHEPROO water-type reader builds it.
        allocate(U_aq(n_p_aq,n_aq_sp))
        U_aq(:,1:n_p_aq)=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq,1:n_p_aq)
        U_aq(:,n_p_aq+1:n_aq_sp)=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq,n_p+1:n_p+n_sec_aq)

        !> Snapshot the current (pre-equilibrium) elemental totals and the conserved components
        allocate(c_aq(n_aq_sp))
        c_aq(1:n_p_aq)=this%concentrations(this%ind_prim_species(1:n_p_aq))
        c_aq(n_p_aq+1:n_aq_sp)=this%concentrations(this%ind_sec_species(1:n_sec_aq))
        T_before=matmul(U_aq,c_aq)
        c_before=this%get_conc_var_act()
        conc_comp=this%compute_conc_comp(c_before) !> u = U*c_nc: invariant under all equilibrium reactions

        !> Re-impose full equilibrium with the reactive zone's minerals/gases (components held fixed)
        allocate(conc_nc(n_v),Delta_n(n_p_aq))
        c1_ig=this%get_c1()
        if (flag_ideal) then
            call this%compute_c_nc_from_u_Newton_ideal(c1_ig,conc_comp,conc_nc,niter,CV_flag)
        else
            c2v_ig=this%get_c2v()
            call this%compute_c_nc_from_u_Newton(c1_ig,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
        end if

        if (.not. CV_flag) then
            Delta_n=0d0 !> equilibrium not reached; report no transfer
            return
        end if

        !> Quantity of each primary species supplied(+)/removed(-) = change in aqueous elemental totals
        c_aq(1:n_p_aq)=this%concentrations(this%ind_prim_species(1:n_p_aq))
        c_aq(n_p_aq+1:n_aq_sp)=this%concentrations(this%ind_sec_species(1:n_sec_aq))
        T_after=matmul(U_aq,c_aq)
        Delta_n=T_after-T_before

        !> Per-phase reaction amounts (mole transfers): compute_Re solves S^T*Re = c_after - c_before
        !> (lambda_r=1) and stores Re on solid_chemistry%rxn_amnt_eq / gas_chemistry%rxn_amnt_eq / this%rxn_amnt_eq (+ rxn_rate_eq).
        if (this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
            call this%compute_Re(c_before,1d0,1d0,0d0)
            !> Update the equilibrium mineral reservoir amounts from re_mean (codebase routine)
            if (this%solid_chemistry%mineral_zone%num_minerals>0) then
                call this%solid_chemistry%compute_conc_minerals_iter(1d0,1d0)
            end if
            !> Update the equilibrium gas reservoir amounts (gas Re and gas concentrations share the
            !> constant-activity-then-variable-activity ordering)
            if (associated(this%gas_chemistry)) then
                do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
                    this%gas_chemistry%concentrations(i)=this%gas_chemistry%concentrations(i)+this%gas_chemistry%rxn_amnt_eq(i)
                    if (this%gas_chemistry%concentrations(i)<0d0) this%gas_chemistry%concentrations(i)=0d0
                end do
            end if
        end if
    end subroutine

    !> \brief Full reaction rates of a past time level, rows ordered as get_S_nc_old:
    !>        [min cst act | gas cst act | aqueous eq | min var act | exch | gas var act | aqueous kin | min kin].
    !> \details Shared body of get_r_old (lev=1: step k), get_r_old_old (lev=2: k-1) and
    !>          get_r_old_old_old (lev=3: k-2). The AQUEOUS reactions (equilibrium and kinetic) belong
    !>          to the OLD solid chemistry, i.e. the reactive zone the water was in at the previous
    !>          step: their count comes from solid_chemistry_old and their rates are carried by the
    !>          water itself (rxn_rate_eq_old*, rxn_rate_kin_old / rk_old_old*). The SOLID and GAS
    !>          reactions (equilibrium minerals, exchange, gases and kinetic minerals) stay in the
    !>          cell: their count AND rates come from the CURRENT solid_chemistry / gas_chemistry.
    !>          The kinetic block is present only when the water has a previous-step kinetic set
    !>          (indices_rk_old set), with the rows of get_rk_old / get_Sk_nc_old.
    function get_r_hist(this,lev) result(r)
        class(aqueous_chemistry_c), intent(in) :: this
        integer(kind=4), intent(in) :: lev !> time level: 1 = k (old), 2 = k-1 (old_old), 3 = k-2 (old_old_old)
        real(kind=8), allocatable :: r(:)
        real(kind=8), allocatable :: re_aq(:) !> aqueous equilibrium rates of that level (carried by the water)
        real(kind=8), allocatable :: re_sol(:) !> solid equilibrium rates of that level (current solid chemistry: [min cst act | min var act | exch])
        real(kind=8), allocatable :: re_gas(:) !> gas equilibrium rates of that level (current gas chemistry: [gas cst act | gas var act])
        real(kind=8), allocatable :: rk_aq(:) !> aqueous kinetic rates of that level (carried by the water)
        real(kind=8), allocatable :: rk_sol(:) !> mineral kinetic rates of that level (current solid chemistry)
        integer(kind=4) :: i,n_r,n_min_cst,n_gas_cst,n_aq_eq,n_min_var,n_exch,n_gas_var,n_aq_kin,n_min_kin

        select case (lev)
        case (1)
            if (allocated(this%rxn_rate_eq_old)) re_aq=this%rxn_rate_eq_old
            if (allocated(this%solid_chemistry%rxn_rate_eq_old)) re_sol=this%solid_chemistry%rxn_rate_eq_old
            if (associated(this%gas_chemistry)) then
                if (allocated(this%gas_chemistry%rxn_rate_eq_old)) re_gas=this%gas_chemistry%rxn_rate_eq_old
            end if
            if (allocated(this%rxn_rate_kin_old)) rk_aq=this%rxn_rate_kin_old
            if (allocated(this%solid_chemistry%rxn_rate_kin_old)) rk_sol=this%solid_chemistry%rxn_rate_kin_old
        case (2)
            if (allocated(this%rxn_rate_eq_old_old)) re_aq=this%rxn_rate_eq_old_old
            if (allocated(this%solid_chemistry%rxn_rate_eq_old_old)) re_sol=this%solid_chemistry%rxn_rate_eq_old_old
            if (associated(this%gas_chemistry)) then
                if (allocated(this%gas_chemistry%rxn_rate_eq_old_old)) re_gas=this%gas_chemistry%rxn_rate_eq_old_old
            end if
            if (allocated(this%rk_old_old)) rk_aq=this%rk_old_old
            if (allocated(this%solid_chemistry%rk_old_old)) rk_sol=this%solid_chemistry%rk_old_old
        case (3)
            if (allocated(this%rxn_rate_eq_old_old_old)) re_aq=this%rxn_rate_eq_old_old_old
            if (allocated(this%solid_chemistry%rxn_rate_eq_old_old_old)) re_sol=this%solid_chemistry%rxn_rate_eq_old_old_old
            if (associated(this%gas_chemistry)) then
                if (allocated(this%gas_chemistry%rxn_rate_eq_old_old_old)) re_gas=this%gas_chemistry%rxn_rate_eq_old_old_old
            end if
            if (allocated(this%rk_old_old_old)) rk_aq=this%rk_old_old_old
            if (allocated(this%solid_chemistry%rk_old_old_old)) rk_sol=this%solid_chemistry%rk_old_old_old
        case default
            error stop "get_r_hist: time level must be 1, 2 or 3"
        end select

        associate(rz => this%solid_chemistry%reactive_zone, rz_old => this%solid_chemistry_old%reactive_zone)
        !> solid & gas reactions: current cell
        n_min_cst=rz%num_minerals_cst_act
        n_gas_cst=rz%gas_phase%num_gases_eq_cst_act
        n_min_var=rz%num_minerals_var_act
        n_exch=rz%cat_exch_zone%num_exch_cats
        n_gas_var=rz%gas_phase%num_gases_eq_var_act
        !> aqueous reactions: old reactive zone
        n_aq_eq=rz_old%chem_syst%num_aq_eq_reacts
        n_aq_kin=0
        n_min_kin=0
        if (this%indices_rk_old%num_cols>0) then !> previous-step kinetic set, as get_rk_old
            n_aq_kin=rz_old%chem_syst%num_aq_kin_reacts
            if (associated(this%solid_chemistry%mineral_zone)) n_min_kin=this%solid_chemistry%mineral_zone%num_minerals_kin
        end if
        allocate(r(n_min_cst+n_gas_cst+n_aq_eq+n_min_var+n_exch+n_gas_var+n_aq_kin+n_min_kin))
        r=0d0
        n_r=0
        do i=1,n_min_cst
            r(n_r+i)=re_sol(i)
        end do
        n_r=n_r+n_min_cst
        do i=1,n_gas_cst
            r(n_r+i)=re_gas(i)
        end do
        n_r=n_r+n_gas_cst
        do i=1,n_aq_eq
            r(n_r+i)=re_aq(i)
        end do
        n_r=n_r+n_aq_eq
        do i=1,n_min_var
            r(n_r+i)=re_sol(n_min_cst+i)
        end do
        n_r=n_r+n_min_var
        do i=1,n_exch
            r(n_r+i)=re_sol(rz%num_minerals+i)
        end do
        n_r=n_r+n_exch
        do i=1,n_gas_var
            r(n_r+i)=re_gas(n_gas_cst+i)
        end do
        n_r=n_r+n_gas_var
        do i=1,n_aq_kin
            r(n_r+i)=rk_aq(i)
        end do
        n_r=n_r+n_aq_kin
        do i=1,n_min_kin
            r(n_r+i)=rk_sol(i)
        end do
        end associate
        !> Faltan gases cineticos
    end function get_r_hist

    !> \brief Full previous-step (k) reaction rates [equilibrium | kinetic], rows ordered as get_S_nc_old.
    !>        Aqueous reactions from the old solid chemistry (rates carried by the water), solid and gas
    !>        reactions from the current solid/gas chemistry; see get_r_hist.
    function get_r_old(this) result(r_old)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old(:)
    r_old=get_r_hist(this,1)
    end function

    !> \brief Returns the OLD reaction rates restricted to AQUEOUS reactions, ordered exactly as the
    !!        rows of get_S_nc_old_aq: [aqueous eq | aqueous kin (linear then redox)].
    !!
    !! \details Aqueous companion of get_r_old for the Lagrangian previous-step split. Sized by
    !! the previous-step reactive zone (solid_chemistry_old), matching get_S_nc_old_aq. The aqueous
    !! equilibrium rates come from this%rxn_rate_eq_old and the aqueous kinetic rates from this%rk_old
    !! (both aqueous-chemistry members, so independent of solid vs solid_old).
    !!
    !! @return r_old_aq aqueous old reaction rates (size = num_aq_eq_reacts + num_aq_kin_reacts)
    !!
    function get_r_old_aq(this) result(r_old_aq)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old_aq(:)
    integer(kind=4) :: i, n_r
    allocate(r_old_aq(this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts+&
        this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts))
    n_r=0
    !> aqueous equilibrium rates
    do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts
        r_old_aq(n_r+i)=this%rxn_rate_eq_old(i)
    end do
    n_r=n_r+this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts
    !> aqueous kinetic rates (linear then redox, as stored in rk_old)
    do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts
        r_old_aq(n_r+i)=this%rxn_rate_kin_old(i)
    end do
    end function

    !> \brief Returns the OLD reaction rates restricted to NON-aqueous reactions, ordered exactly as
    !!        the rows of get_S_nc_non_aq:
    !!        [min cst act | gas eq cst act | min eq var act | exch | gas eq var act | mineral kin].
    !!
    !! \details Non-aqueous companion of get_r_old for the Lagrangian previous-step split. It is
    !! evaluated at the CURRENT cell (uses this%solid_chemistry / this%gas_chemistry), matching
    !! get_S_nc_non_aq, because the solid & gas reactions are "fixed in space" at the current cell. The
    !! VALUES are the previous-step (k) rates: rxn_rate_eq_old for the equilibrium reactions and rk_old for
    !! the mineral kinetic reactions. This is the same source get_r_old uses for its solid/gas rows
    !! (see get_r_hist), so both agree row by row on the non-aqueous reactions.
    !!
    !! @return r_old_non_aq non-aqueous old reaction rates (current reactive zone)
    !!
    function get_r_old_non_aq(this) result(r_old_non_aq)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old_non_aq(:)
    integer(kind=4) :: i, n_r
    allocate(r_old_non_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts+&
        this%solid_chemistry%mineral_zone%num_minerals_kin))
    n_r=0
    !> constant activity minerals (equilibrium)
    do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
        r_old_non_aq(n_r+i)=this%solid_chemistry%rxn_rate_eq_old(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%num_minerals_cst_act
    !> constant activity gases (equilibrium)
    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        r_old_non_aq(n_r+i)=this%gas_chemistry%rxn_rate_eq_old(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
    !> variable activity minerals (equilibrium)
    do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
        r_old_non_aq(n_r+i)=this%solid_chemistry%rxn_rate_eq_old(&
            this%solid_chemistry%reactive_zone%num_minerals_cst_act+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%num_minerals_var_act
    !> exchange (surface complexes)
    do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
        r_old_non_aq(n_r+i)=this%solid_chemistry%rxn_rate_eq_old(&
            this%solid_chemistry%reactive_zone%num_minerals+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
    !> variable activity gases (equilibrium)
    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
        r_old_non_aq(n_r+i)=this%gas_chemistry%rxn_rate_eq_old(&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
    !> mineral kinetic
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        r_old_non_aq(n_r+i)=this%solid_chemistry%rxn_rate_kin_old(i)
    end do
    end function

    !> \brief Returns the FULL reaction rates from two time steps ago, ordered exactly as get_r_old.
    !!
    !! \details Same aqueous/solid split as get_r_old (see get_r_hist) with the two-steps-ago history:
    !! rxn_rate_eq_old_old for the equilibrium reactions and rk_old_old for the kinetic reactions.
    !! Used by the downstream first-derivative extrapolation (option 2), which estimates
    !! r^{k+1} = (1+mu) r^k - mu r^{k-1} for BOTH equilibrium and kinetic reactions.
    function get_r_old_old(this) result(r_old_old)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old_old(:)
    r_old_old=get_r_hist(this,2)
    end function

    !> \brief Returns the FULL reaction rates from three time steps ago, ordered exactly as get_r_old.
    !!
    !! \details Same aqueous/solid split as get_r_old (see get_r_hist) with the three-steps-ago history:
    !! rxn_rate_eq_old_old_old for the equilibrium reactions and rk_old_old_old for the kinetic reactions.
    !! Used by the downstream second-derivative extrapolation (option 3), which estimates
    !! r^{k+1} = 3 r^k - 3 r^{k-1} + r^{k-2}.
    function get_r_old_old_old(this) result(r_old_old_old)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old_old_old(:)
    r_old_old_old=get_r_hist(this,3)
    end function
    
    function get_rxn_amnt_kin(this) result(Rk)
    !!> This function returns the kinetic reaction amounts ordered in aqueous & mineral reactions
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: Rk(:)
    integer(kind=4) :: i
    allocate(Rk(this%indices_rk%num_cols))
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        Rk(i)=this%rxn_amnt_kin(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        Rk(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rxn_amnt_kin(i)
    end do
    !> Faltan gases
    end function

    !> \brief Equilibrium reaction amounts of the step, ordered as the rows of the reactive zone's
    !!        stoichiometric matrix (the Se_nc that compute_contrib_rxn_amnt_eq pairs them with).
    !! \details Inverse of the distribution done by compute_Re: constant-activity minerals, constant-activity
    !!          gases, aqueous reactions, variable-activity minerals, exchange cations, variable-activity
    !!          gases. The gas entries stay 0 when the water has no gas chemistry.
    function get_rxn_amnt_eq(this) result(Re)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: Re(:)
    integer(kind=4) :: i,n
    associate(rz => this%solid_chemistry%reactive_zone)
    allocate(Re(rz%speciation_alg%num_eq_reactions))
    Re=0d0
    n=0
    do i=1,rz%num_minerals_cst_act
        Re(n+i)=this%solid_chemistry%rxn_amnt_eq(i)
    end do
    n=n+rz%num_minerals_cst_act
    if (associated(this%gas_chemistry)) then
        do i=1,rz%gas_phase%num_gases_eq_cst_act
            Re(n+i)=this%gas_chemistry%rxn_amnt_eq(i)
        end do
    end if
    n=n+rz%gas_phase%num_gases_eq_cst_act
    do i=1,rz%chem_syst%num_aq_eq_reacts
        Re(n+i)=this%rxn_amnt_eq(i)
    end do
    n=n+rz%chem_syst%num_aq_eq_reacts
    do i=1,rz%num_minerals_var_act
        Re(n+i)=this%solid_chemistry%rxn_amnt_eq(rz%num_minerals_cst_act+i)
    end do
    n=n+rz%num_minerals_var_act
    do i=1,rz%cat_exch_zone%num_exch_cats
        Re(n+i)=this%solid_chemistry%rxn_amnt_eq(rz%num_minerals+i)
    end do
    n=n+rz%cat_exch_zone%num_exch_cats
    if (associated(this%gas_chemistry)) then
        do i=1,rz%gas_phase%num_gases_eq_var_act
            Re(n+i)=this%gas_chemistry%rxn_amnt_eq(rz%gas_phase%num_gases_eq_cst_act+i)
        end do
    end if
    end associate
    end function

    !> \brief Previous-step kinetic reaction rates of this water in its OLD kinetic reaction set.
    !> \details The previous-step rates of the AQUEOUS kinetic reactions belong to the old solid
    !>          chemistry (the reactive zone the water was in at the previous step) and are carried
    !>          by the water itself (rxn_rate_kin_old). The previous-step rates of the SOLID kinetic
    !>          reactions (minerals) belong to the CURRENT solid chemistry: the minerals stay in the
    !>          cell, so the old set and the new set (get_rk_old_new_set) share them
    !>          (solid_chemistry%rxn_rate_kin_old) and differ only in their aqueous part.
    !>          Rows ordered as get_Sk_nc_old (aqueous kinetic reactions of the old reactive zone,
    !>          then the kinetic minerals of the current solid chemistry), the stoichiometry this
    !>          vector must be paired with (previous-step kinetic term of the mixers). Zero-sized
    !>          when the water has no previous-step kinetic set (indices_rk_old unset), which then
    !>          contributes no previous-step kinetic term.
    function get_rk_old(this) result(rk_old)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: rk_old(:)
        integer(kind=4) :: i,n_aq_kin,n_min_kin
        if (this%indices_rk_old%num_cols==0) then
            allocate(rk_old(0))
            return
        end if
        n_aq_kin=this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts
        n_min_kin=0
        if (associated(this%solid_chemistry%mineral_zone)) n_min_kin=this%solid_chemistry%mineral_zone%num_minerals_kin
        if (this%indices_rk_old%num_cols/=n_aq_kin+n_min_kin) then
            print '(A,I0,A,I0,A,I0,A,I0)', "ERROR get_rk_old: water ", this%id, ": indices_rk_old%num_cols=", this%indices_rk_old%num_cols, &
                " /= aqueous kinetic of the old reactive zone (", n_aq_kin, ") + kinetic minerals of the current solid chemistry (", &
                n_min_kin
            error stop "get_rk_old: old kinetic reaction set inconsistent with the old reactive zone / current solid chemistry"
        end if
        allocate(rk_old(n_aq_kin+n_min_kin))
        rk_old=0d0
        do i=1,n_aq_kin
            rk_old(i)=this%rxn_rate_kin_old(i)
        end do
        do i=1,n_min_kin
            rk_old(n_aq_kin+i)=this%solid_chemistry%rxn_rate_kin_old(i)
        end do
        !> Faltan gases
    end function

    !> \brief Previous-step kinetic reaction rates of this water in its NEW (current) kinetic
    !>        reaction set. As in get_rk_old, the old rates of the AQUEOUS kinetic reactions belong to
    !>        the OLD solid chemistry (the reactive zone the water was in at the previous step; rates
    !>        carried by the water, rxn_rate_kin_old), and those of the SOLID kinetic reactions
    !>        (minerals) to the CURRENT solid chemistry (solid_chemistry%rxn_rate_kin_old: the minerals
    !>        stay in the cell, so both sets share them). The two getters differ in the kinetic set the
    !>        vector is validated against: here the new set (indices_rk), the rows of get_Sk_nc / rk_new.
    !> \details This is the vector combined with the current-set stoichiometry and rates: the Newton
    !>          predictor and the mean rate rk_avg = theta*rk_new + (1-theta)*rk_old of the current
    !>          cell's minerals, the lumped r_tilde and the Euler explicit reaction term (U_SkT_prod
    !>          of the current reactive zone). Zero-sized when the water has no kinetic reactions
    !>          (indices_rk unset).
    function get_rk_old_new_set(this) result(rk_old)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: rk_old(:)
        integer(kind=4) :: i,n_aq_kin,n_min_kin
        if (this%indices_rk%num_cols==0) then
            allocate(rk_old(0))
            return
        end if
        n_aq_kin=this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts !> aqueous kinetic reactions: old solid chemistry (old reactive zone)
        n_min_kin=0
        if (associated(this%solid_chemistry%mineral_zone)) n_min_kin=this%solid_chemistry%mineral_zone%num_minerals_kin !> solid kinetic reactions: current solid chemistry
        if (this%indices_rk%num_cols/=n_aq_kin+n_min_kin) then
            print '(A,I0,A,I0,A,I0,A,I0)', "ERROR get_rk_old_new_set: water ", this%id, ": indices_rk%num_cols=", &
                this%indices_rk%num_cols, " /= aqueous kinetic of the old reactive zone (", n_aq_kin, &
                ") + kinetic minerals of the current solid chemistry (", n_min_kin
            error stop "get_rk_old_new_set: new kinetic reaction set inconsistent with the old reactive zone / current solid chemistry"
        end if
        allocate(rk_old(n_aq_kin+n_min_kin))
        rk_old=0d0
        do i=1,n_aq_kin
            rk_old(i)=this%rxn_rate_kin_old(i) !> aqueous: rates carried by the water
        end do
        do i=1,n_min_kin
            rk_old(n_aq_kin+i)=this%solid_chemistry%rxn_rate_kin_old(i) !> minerals: rates of the current cell's solid chemistry
        end do
        !> Faltan gases
    end function get_rk_old_new_set

    !> \brief Returns the variable activity stoichiometric matrix of kinetic reactions
    !!
    !! @return Sk_nc kinetic stoichiometric matrix
    !!
    function get_Sk_nc(this) result(Sk_nc) !> returns variable activity stoichiometric matrix of kinetic reactions
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: Sk_nc(:,:)
        integer(kind=4) :: i,num_rk
        allocate(Sk_nc(this%indices_rk%num_cols,this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species))
        num_rk=0
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            Sk_nc(num_rk+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            Sk_nc(num_rk+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            Sk_nc(num_rk+i,:)=&
                this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        !> Faltan gases
    end function

    !> \brief Returns the variable activity stoichiometric matrix for all reactions
    !!
    !! \details Assembles the full stoichiometric matrix combining equilibrium
    !! and kinetic reactions for variable activity species.
    !!
    !! @param ind_var_act_species optional indices of variable activity species
    !! @return S_nc full stoichiometric matrix
    !!
    function get_S_nc(this,ind_var_act_species) result(S_nc) !> returns variable activity stoichiometric matrix of kinetic reactions
        class(aqueous_chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_var_act_species(:) !< Indices of variable activity species in chemical system (required: the equilibrium block is built from these columns)
        real(kind=8), allocatable :: S_nc(:,:)
        integer(kind=4) :: i,j,n_r,num_rk
        allocate(S_nc(this%indices_rk%num_cols+this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species))
        ! S_nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,:)=&
        !     this%solid_chemistry%reactive_zone%get_Se_nc_react_zone(ind_var_act_species)
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions     
            do j=1,size(ind_var_act_species)
                !< Extract column for each specified species
                S_nc(i,j)=this%solid_chemistry%reactive_zone%stoich_mat(&
                    i,this%solid_chemistry%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        num_rk=0
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            S_nc(n_r+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        n_r=n_r+num_rk
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            S_nc(n_r+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        n_r=n_r+num_rk
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            S_nc(n_r+i,:)=&
                this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        !> Faltan gases
    end function

    !> \brief Returns the FULL variable activity stoichiometric matrix of the previous-step reaction set,
    !!        rows ordered exactly as get_r_old (get_r_hist):
    !!        [min cst act | gas cst act | aqueous eq | min var act | exch | gas var act | aqueous kin | min kin].
    !!
    !! \details Same aqueous/solid split as get_r_old: the AQUEOUS equilibrium rows come from the OLD
    !! reactive zone (solid_chemistry_old, the zone the water was in at the previous step) and the
    !! SOLID and GAS equilibrium rows (constant- and variable-activity minerals, exchange, gases) from
    !! the CURRENT reactive zone, since those phases stay in the cell. The kinetic rows are
    !! get_Sk_nc_old (aqueous kinetic reactions of the old zone, then the kinetic minerals of the
    !! current solid chemistry; absent when the water has no previous-step kinetic set). Columns are
    !! the variable-activity species of the current reactive zone, the basis of the mixed
    !! concentrations the previous-step term is added to.
    !!
    !! @param ind_var_act_species indices of variable activity species
    !! @return S_nc_old full stoichiometric matrix of the previous-step reaction set
    !!
    function get_S_nc_old(this,ind_var_act_species) result(S_nc_old)
        class(aqueous_chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_var_act_species(:) !< Indices of variable activity species in chemical system (required: the equilibrium block is built from these columns)
        real(kind=8), allocatable :: S_nc_old(:,:)
        real(kind=8), allocatable :: Sk_nc_old(:,:) !> kinetic rows (get_Sk_nc_old)
        integer(kind=4) :: i,j,n_r,n_v,n_cols,n_cst,n_aq_eq,n_var_non_aq,off_old,off_cur
        associate(rz => this%solid_chemistry%reactive_zone, rz_old => this%solid_chemistry_old%reactive_zone)
        n_v=rz%speciation_alg%num_var_act_species
        n_cols=min(n_v,size(ind_var_act_species))
        n_cst=rz%num_minerals_cst_act+rz%gas_phase%num_gases_eq_cst_act !> constant-activity mineral & gas rows (current zone)
        n_aq_eq=rz_old%chem_syst%num_aq_eq_reacts !> aqueous equilibrium rows (old zone)
        n_var_non_aq=rz%num_minerals_var_act+rz%cat_exch_zone%num_exch_cats+rz%gas_phase%num_gases_eq_var_act !> variable-activity mineral, exchange & gas rows (current zone)
        Sk_nc_old=this%get_Sk_nc_old()
        allocate(S_nc_old(n_cst+n_aq_eq+n_var_non_aq+size(Sk_nc_old,1),n_v))
        S_nc_old=0d0
        !> rows of a reactive zone's stoich_mat are its ind_eq_reacts order:
        !> [cst act minerals & gases | aqueous eq | var act minerals, exchange & gases]
        n_r=0
        do i=1,n_cst
            do j=1,n_cols
                S_nc_old(n_r+i,j)=rz%stoich_mat(i,rz%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=n_r+n_cst
        off_old=rz_old%num_minerals_cst_act+rz_old%gas_phase%num_gases_eq_cst_act
        do i=1,n_aq_eq
            do j=1,n_cols
                S_nc_old(n_r+i,j)=rz_old%stoich_mat(off_old+i,rz_old%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=n_r+n_aq_eq
        off_cur=n_cst+rz%chem_syst%num_aq_eq_reacts
        do i=1,n_var_non_aq
            do j=1,n_cols
                S_nc_old(n_r+i,j)=rz%stoich_mat(off_cur+i,rz%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=n_r+n_var_non_aq
        if (size(Sk_nc_old,1)>0) S_nc_old(n_r+1:,:)=Sk_nc_old !> same columns: get_Sk_nc_old uses the current zone's variable-activity species too
        end associate
    end function

    !> \brief Returns the variable activity stoichiometric matrix from the previous time step restricted to aqueous reactions
    !!
    !! @param ind_var_act_species optional indices of variable activity species
    !! @return S_nc_old stoichiometric matrix from previous time step
    !!
    function get_S_nc_old_aq(this,ind_var_act_species) result(S_nc_old_aq) !> returns variable activity stoichiometric matrix of kinetic reactions
        class(aqueous_chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_var_act_species(:) !< Indices of variable activity species in chemical system (required: the equilibrium block is built from these columns)
        real(kind=8), allocatable :: S_nc_old_aq(:,:)
        integer(kind=4) :: i,n_r,num_rk, j
        allocate(S_nc_old_aq(this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts+&
            this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts,&
            this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species))
        ! S_nc_old(1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions,:)=&
        !     this%solid_chemistry_old%reactive_zone%get_Se_nc_react_zone(ind_var_act_species)
        ! do i=1,this%solid_chemistry%reactive_zone%num_non_flow_species     
        !     do j=1,size(ind_var_act_species)
        !         !< Extract column for each specified species
        !         S_nc_old_aq(i,j)=this%solid_chemistry%reactive_zone%stoich_mat(&
        !             i,this%solid_chemistry%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
        !     end do
        ! end do
        ! n_r=this%solid_chemistry%reactive_zone%num_non_flow_species
        !> AQUEOUS EQUILIBRIUM block. In the reactive-zone stoich_mat the equilibrium reactions are
        !> ordered [min cst act | gas eq cst act | AQUEOUS eq | min var act | exch | gas eq var act]
        !> (same order as get_r_old), so the aqueous eq rows are NOT the first rows: they start after
        !> the constant-activity mineral and gas equilibrium reactions.
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts     
            do j=1,size(ind_var_act_species)
                !< Extract column for each specified species
                S_nc_old_aq(i,j)=this%solid_chemistry_old%reactive_zone%stoich_mat(&
                    this%solid_chemistry_old%reactive_zone%num_minerals_cst_act+&
                    this%solid_chemistry_old%reactive_zone%gas_phase%num_gases_eq_cst_act+i,&
                    this%solid_chemistry_old%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts   
        ! num_rk=0
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
            S_nc_old_aq(n_r+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        end do
        n_r=n_r+this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
        !> Redox-kinetic rows follow the linear-kinetic rows in chem_syst%Sk, so offset by num_lin_kin_reacts.
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
            S_nc_old_aq(n_r+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(&
                this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts+i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        end do
        n_r=n_r+this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
        ! do i=1,this%solid_chemistry_old%mineral_zone%num_minerals_kin
        !     S_nc_old_aq(n_r+i,:)=&
        !         this%solid_chemistry_old%reactive_zone%chem_syst%Sk(i,&
        !         this%solid_chemistry_old%mineral_zone%ind_min_chem_syst(i),&
        !         1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        ! end do
    end function

    !> \brief Returns the variable activity stoichiometric matrix (CURRENT/new reactive zone)
    !!        restricted to NON-aqueous reactions (solid & gas equilibrium + mineral kinetic).
    !!
    !! \details Companion of get_S_nc_old_aq for the Lagrangian previous-step split: aqueous
    !! reactions are taken from where the water WAS at the previous step (get_S_nc_old_aq, which
    !! uses solid_chemistry_old), whereas solid & gas reactions are taken from where the water IS
    !! at the current step -- so this function uses the CURRENT solid_chemistry (the "new"
    !! stoichiometric matrix), NOT solid_chemistry_old.
    !!
    !! Rows are the NON-aqueous subset of get_S_nc, in the same order:
    !!   [min cst act | gas eq cst act | min eq var act | exch | gas eq var act | mineral kin]
    !! i.e. every equilibrium reaction EXCEPT the aqueous eq block (which sits between the leading
    !! cst-act min/gas rows and the var-act min/exch/gas rows in stoich_mat), followed by the
    !! mineral kinetic reactions (the aqueous kinetic block -- linear + redox -- is skipped in Sk).
    !!
    !! @param ind_var_act_species indices of variable activity species in chemical system
    !! @return S_nc_non_aq non-aqueous stoichiometric matrix from the current time step
    !!
    function get_S_nc_non_aq(this,ind_var_act_species) result(S_nc_non_aq) !> returns non-aqueous variable activity stoichiometric matrix (current reactive zone)
        class(aqueous_chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_var_act_species(:) !< Indices of variable activity species in chemical system (the equilibrium block is built from these columns)
        real(kind=8), allocatable :: S_nc_non_aq(:,:)
        integer(kind=4) :: i,n_r,j,n_head,n_tail
        !> n_head = leading non-aqueous equilibrium reactions (cst-act minerals + cst-act gases),
        !> which precede the aqueous eq block in stoich_mat.
        n_head=this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        !> n_tail = trailing non-aqueous equilibrium reactions (var-act minerals + exchange +
        !> var-act gases), which follow the aqueous eq block in stoich_mat.
        n_tail=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-n_head-&
            this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
        allocate(S_nc_non_aq(n_head+n_tail+this%solid_chemistry%mineral_zone%num_minerals_kin,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
        n_r=0
        !> NON-AQUEOUS EQUILIBRIUM block 1: leading cst-act minerals & gases (stoich_mat rows 1..n_head)
        do i=1,n_head
            do j=1,size(ind_var_act_species)
                !< Extract column for each specified species
                S_nc_non_aq(n_r+i,j)=this%solid_chemistry%reactive_zone%stoich_mat(&
                    i,this%solid_chemistry%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=n_r+n_head
        !> NON-AQUEOUS EQUILIBRIUM block 2: var-act minerals, exchange, var-act gases (stoich_mat rows
        !> after the aqueous eq block: n_head+num_aq_eq_reacts+1 .. num_eq_reactions).
        do i=1,n_tail
            do j=1,size(ind_var_act_species)
                !< Extract column for each specified species
                S_nc_non_aq(n_r+i,j)=this%solid_chemistry%reactive_zone%stoich_mat(&
                    n_head+this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts+i,&
                    this%solid_chemistry%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=n_r+n_tail
        !> NON-AQUEOUS KINETIC block: mineral kinetic reactions (chem_syst%Sk rows after the aqueous
        !> kinetic block -- linear + redox = num_aq_kin_reacts -- selected by ind_min_chem_syst).
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            S_nc_non_aq(n_r+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(&
                this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)
        end do
    end function

    !> \brief Returns the kinetic stoichiometric matrix of the OLD kinetic reaction set
    !!
    !! \details Rows ordered as get_rk_old: the aqueous kinetic reactions (linear, redox) of the OLD
    !!          reactive zone (solid_chemistry_old), then the kinetic minerals of the CURRENT solid
    !!          chemistry (they stay in the cell, so the old set shares them with the new set).
    !!          Columns are the variable-activity species of the current reactive zone, the basis of
    !!          the mixed concentrations the previous-step kinetic term is added to.
    !! @return Sk_nc_old kinetic stoichiometric matrix of the old kinetic set
    !!
    function get_Sk_nc_old(this) result(Sk_nc_old)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: Sk_nc_old(:,:)
        integer(kind=4) :: i,num_rk,n_lin,n_redox,n_min_kin,n_v
        n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        !> no previous-step kinetic set: no rows, as get_rk_old
        if (this%indices_rk_old%num_cols==0) then
            allocate(Sk_nc_old(0,n_v))
            return
        end if
        n_lin=this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
        n_redox=this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
        n_min_kin=0
        if (associated(this%solid_chemistry%mineral_zone)) n_min_kin=this%solid_chemistry%mineral_zone%num_minerals_kin
        !> rows = the old kinetic set; a mismatch with indices_rk_old would leave unset (garbage)
        !> rows, so it is an error
        if (this%indices_rk_old%num_cols/=n_lin+n_redox+n_min_kin) then
            print '(A,I0,A,I0,A,I0,A,I0)', "ERROR get_Sk_nc_old: water ", this%id, ": indices_rk_old%num_cols=", &
                this%indices_rk_old%num_cols, " /= aqueous kinetic of the old reactive zone (", n_lin+n_redox, &
                ") + kinetic minerals of the current solid chemistry (", n_min_kin
            error stop "get_Sk_nc_old: old kinetic reaction set inconsistent with the old reactive zone / current solid chemistry"
        end if
        allocate(Sk_nc_old(n_lin+n_redox+n_min_kin,n_v))
        Sk_nc_old=0d0
        num_rk=0
        do i=1,n_lin
            Sk_nc_old(num_rk+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,1:n_v)
        end do
        num_rk=num_rk+n_lin
        do i=1,n_redox
            Sk_nc_old(num_rk+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,1:n_v)
        end do
        num_rk=num_rk+n_redox
        do i=1,n_min_kin
            Sk_nc_old(num_rk+i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(num_rk+&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i),1:n_v)
        end do
        !> Faltan gases
    end function
    
    !> \brief Computes kinetic reaction amounts during a time step
    !!
    !! \details First computes the time-weighted average rate rk_avg = theta*rk_new + (1-theta)*rk_old
    !!          (stored in rxn_rate_kin_avg), then the reaction amounts Rk = rk_avg*Delta_t.
    !!
    !! @param theta time weighting factor
    !! @param Delta_t time step
    !! @param lambda_r reaction mixing ratio
    !!
    subroutine compute_rxn_amnt_kin(this,Delta_t) !> computes kinetic reaction amounts during a time step
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    !real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio
    integer(kind=4) :: i
    ! !> Time-weighted average kinetic reaction rate: rk_avg = theta*rk_new + (1-theta)*rk_old
    ! do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
    !     this%rxn_rate_kin_avg(i)=theta*this%rxn_rate_kin(i)+(1d0-theta)*this%rxn_rate_kin_old(i)
    ! end do
    ! do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
    !     this%solid_chemistry%rxn_rate_kin_avg(i)=theta*this%solid_chemistry%rxn_rate_kin(i)+(1d0-theta)*this%solid_chemistry%rxn_rate_kin_old(i)
    ! end do
    ! do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    !     this%rxn_rate_kin_avg(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
    !         theta*this%rxn_rate_kin(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)+&
    !         (1d0-theta)*this%rxn_rate_kin_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
    ! end do
    !> Kinetic reaction amounts over the step: Rk = rk_avg*Delta_t
    this%rxn_amnt_kin=this%rxn_rate_kin_avg*Delta_t
    this%solid_chemistry%rxn_amnt_kin=this%solid_chemistry%rxn_rate_kin_avg*Delta_t
    !> Contribution of kinetic reaction amounts to species mass balance: rho_kin = lambda_r*Sk_nc^T*Rk
    ! if (this%indices_rk%num_cols>0) then
    !     this%contrib_rxn_amnt_kin=lambda_r*matmul(transpose(this%get_Sk_nc()),this%get_Rk())
    ! else
    !     this%contrib_rxn_amnt_kin=[(0d0,i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)]
    ! end if
    end subroutine

    !> \brief Computes kinetic reaction amounts during a time step
    !!
    !! \details First computes the time-weighted average rate rk_avg = theta*rk_new + (1-theta)*rk_old
    !!          (stored in rxn_rate_kin_avg), then the reaction amounts Rk = rk_avg*Delta_t.
    !!
    !! @param theta time weighting factor
    !! @param Delta_t time step
    !! @param lambda_r reaction mixing ratio
    !!
    subroutine compute_contrib_rxn_amnt_kin(this,lambda_r) !> computes kinetic reaction amounts during a time step
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: theta !> time weighting factor
    !real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio
    integer(kind=4) :: i
    ! !> Time-weighted average kinetic reaction rate: rk_avg = theta*rk_new + (1-theta)*rk_old
    ! do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
    !     this%rxn_rate_kin_avg(i)=theta*this%rxn_rate_kin(i)+(1d0-theta)*this%rxn_rate_kin_old(i)
    ! end do
    ! do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
    !     this%solid_chemistry%rxn_rate_kin_avg(i)=theta*this%solid_chemistry%rxn_rate_kin(i)+(1d0-theta)*this%solid_chemistry%rxn_rate_kin_old(i)
    ! end do
    ! do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    !     this%rxn_rate_kin_avg(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
    !         theta*this%rxn_rate_kin(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)+&
    !         (1d0-theta)*this%rxn_rate_kin_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
    ! end do
    ! !> Kinetic reaction amounts over the step: Rk = rk_avg*Delta_t
    ! this%rxn_amnt_kin=this%rxn_rate_kin_avg*Delta_t
    ! this%solid_chemistry%rxn_amnt_kin=this%solid_chemistry%rxn_rate_kin_avg*Delta_t
    !> Contribution of kinetic reaction amounts to species mass balance: rho_kin = lambda_r*Sk_nc^T*Rk
    ! if (this%indices_rk%num_cols>0) then
    !     this%contrib_rxn_amnt_kin=lambda_r*matmul(transpose(this%get_Sk_nc()),this%get_Rk())
    ! else
    !     this%contrib_rxn_amnt_kin=[(0d0,i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)]
    ! end if
    this%contrib_rxn_amnt_kin=lambda_r*matmul(transpose(this%get_Sk_nc()),this%get_rxn_amnt_kin())
    end subroutine

    !> \brief Computes the contribution of the equilibrium reaction amounts to the species mass balance
    !!
    !! \details rho_eq = lambda_r*Se_nc^T*Re, the equilibrium counterpart of compute_contrib_rxn_amnt_kin:
    !!          Se_nc is the variable-activity block of the reactive zone's equilibrium stoichiometric
    !!          matrix (stoich_mat(:,1:n_v), referenced in place as in compute_Re) and Re the equilibrium
    !!          reaction amounts of the step in the same row order (get_rxn_amnt_eq). Stored in
    !!          contrib_rxn_amnt_eq (size num_var_act_species). Note that compute_Re builds rho_eq the
    !!          other way round (from the speciated and mixed concentrations) and then solves for Re; this
    !!          routine recomputes it from the stored Re (e.g. after the amounts have been modified).
    !!
    !! @param lambda_r reaction mixing ratio
    !!
    subroutine compute_contrib_rxn_amnt_eq(this,lambda_r) !> computes equilibrium reaction amount contributions to the species mass balance
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio
    !> rho_eq = lambda_r*Se_nc^T*Re; matmul(vector, matrix) = Se_nc^T*Re without forming the transpose
    this%contrib_rxn_amnt_eq=lambda_r*matmul(this%get_rxn_amnt_eq(),this%solid_chemistry%reactive_zone%stoich_mat(:,&
        1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    end subroutine

    !> \brief Sets the mean kinetic reaction rates during a time step
    !!
    !! @param rk_mean mean kinetic reaction rates during a time step
    !!
    subroutine set_rxn_rate_kin_avg(this,rk_mean)
        !!> This subroutine sets the mean kinetic reaction rates during a time step
        !! We assume that rk_mean is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: rk_mean(:) !> mean kinetic reaction rates during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%rxn_rate_kin_avg(i)=rk_mean(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rxn_rate_kin_avg(i)=rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%rxn_rate_kin_avg(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
        !this%Rk_mean=this%Rk_mean*Delta_t
        !this%solid_chemistry%Rk_mean=this%solid_chemistry%Rk_mean*Delta_t
    end subroutine

    !> \brief Sets the kinetic reaction amounts during a time step
    !!
    !! @param Rk kinetic reaction amounts during a time step
    !!
    subroutine set_Rk(this,Rk)
        !!> This subroutine sets the kinetic reaction amounts during a time step
        !! We assume that Rk is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: Rk(:) !> kinetic reaction amounts during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%rxn_amnt_kin(i)=Rk(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rxn_amnt_kin(i)=Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%rxn_amnt_kin(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
        !this%Rk_mean=this%Rk_mean*Delta_t
        !this%solid_chemistry%Rk_mean=this%solid_chemistry%Rk_mean*Delta_t
    end subroutine

    !> \brief Sets the old kinetic reaction rates
    !!
    !! @param rk_old old kinetic reaction rates
    !!
    subroutine set_rk_old(this,rk_old)
        !!> This subroutine sets the old kinetic reaction rates (chapuza) during a time step
        !! We assume that rk_old is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: rk_old(:) !> old kinetic reaction rates during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
            this%rxn_rate_kin_old(i)=rk_old(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rxn_rate_kin_old(i)=rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)
        end do
        !> Faltan gases
    end subroutine
    
    !> \brief Returns the aqueous component concentrations
    !!
    !! @return u_aq aqueous component concentrations
    !!
    function get_u_aq(this) result(u_aq)
    !> This function returns the aqueous component concentrations
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: u_aq(:)
    allocate(u_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    u_aq=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq,&
        this%concentrations(this%ind_var_act_species))
    end function get_u_aq

    !> \brief Sets speciation algebra variables in the reactive zone
    !!
    !! \details Configures the speciation algebra object, including component
    !! matrix definition, species indices, stoichiometric matrix, and
    !! equilibrium reactions.
    !!
    !! @param flag_comp flag indicating component matrix definition
    !! @param flag_surf output flag for surface complexation reactions
    !! @param flag_Se output flag indicating if stoichiometric matrix was modified
    !! @param swap output indices to swap in stoichiometric matrix
    !! @param ind_eq_reacts optional indices of equilibrium reactions
    !!
    subroutine set_spec_alg_aq_chem(this,flag_comp,flag_surf,flag_Se,swap,ind_eq_reacts)
        !!> This subroutine sets all necessary variables in the speciation algebra object of the reactive zone pointer
        class(aqueous_chemistry_c) :: this
        logical, intent(in) :: flag_comp !> flag to indicate component matrix definition
        logical, intent(out) :: flag_surf !> flag to indicate if there are surface complexation reactions
        logical, intent(out) :: flag_Se !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
        integer(kind=4), intent(out) :: swap(:) !> species or reactions to swap in stoichiometric matrix (already allocated)
        integer(kind=4), intent(in), optional :: ind_eq_reacts(:) !> indices of equilibrium reactions in chemical system
        !logical :: flag_surf !> flag to indicate if there are surface complexation reactions

        call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(flag_comp)
        write(0,'(A,L1)') ' [SSA-DBG] enter set_spec_alg assoc_sc=', associated(this%solid_chemistry)
        if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
            flag_surf=.true.
        else
            flag_surf=.false.
        end if
        !> Ensure solid-chemistry arrays are allocated for BOTH cases: the no-exchange case
        !> and the exchange case (get_c1's surface-site term and compute_conc_surf_compl write
        !> into solid_chemistry%concentrations, so it must be sized num_solids beforehand).
        call this%solid_chemistry%allocate_conc_solids()
        write(0,'(A,L1,A,I0)') ' [SSA-DBG] post allocate_conc_solids alloc=', &
            allocated(this%solid_chemistry%concentrations), ' nsolids=', this%solid_chemistry%num_solids
        call this%solid_chemistry%allocate_activities()
        call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
        call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
        call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
        call this%set_ind_species()
        call this%solid_chemistry%reactive_zone%set_ind_eq_reacts(ind_eq_reacts)
        call this%solid_chemistry%reactive_zone%allocate_ind_var_act_species()
        call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
        call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
        call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
        !allocate(swap(2))
        !> compute_speciation_alg_arrays builds ONLY the structural, K-independent arrays
        !> on the (possibly shared) reactive zone. The aqueous chemistry owns both the solid
        !> (reactive zone) and gas chemistry, so it then builds the modified equilibrium
        !> constants K' = K/(prod a_cst) (constant-activity gas partial pressures fold in;
        !> minerals have a=1) and its transformed constants logK_tilde/logK_star PER-WATER.
        block
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
            !> Build this water's per-reaction T-dependent logK FIRST (compute_logK_star consumes it
            !> via get_eq_csts_aq_chem = 10**logK), then its transformed constants logK_tilde/logK_star
            !> (which also fold in this water's constant-activity gas partial pressures).
            call this%compute_logK_PHREEQC_aq_chem() !> build this water's per-reaction logK at its own temperature
            if (flag_comp) then
                call this%compute_logK_star() !> component elimination (flag_comp=.true.): build logK_star
            else
                call this%compute_logK_tilde_aq_chem() !> all secondary species (flag_comp=.false.): build logK_tilde
            end if
        end block
    end subroutine

    !> \brief Eliminates constant activity species from the component matrix
    !!
    !! \details Reconfigures the speciation algebra by removing constant
    !! activity species and recomputing the component matrix and stoichiometric
    !! arrays.
    !!
    subroutine change_spec_alg_aq_chem(this)
    !!> This subroutine eliminates the constant activity species from the component matrix in the speciation algebra object of the &
        !! reactive zone pointer
        class(aqueous_chemistry_c) :: this !> aqueous chemistry object

        logical :: flag_Se !> flag to indicate if stoichiometric matrix has been modified
        integer(kind=4), allocatable :: swap(:) !> species or reactions to swap in stoichiometric matrix

        call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(.true.)
        call this%set_ind_species()
        !> The chemical system was just rearranged by change_spec_alg_chem_syst,
        !> so the previously stored ind_non_flow_species point to stale chem_syst
        !> positions. Regenerate them from the (now updated) chem_syst speciation
        !> layout before using them in set_ind_eq_reacts.
        call this%solid_chemistry%reactive_zone%allocate_ind_non_flow_species()
        call this%solid_chemistry%reactive_zone%set_ind_non_flow_species()
        call this%solid_chemistry%reactive_zone%set_ind_eq_reacts()
        call this%solid_chemistry%reactive_zone%allocate_ind_var_act_species()
        call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
        call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
        call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
        allocate(swap(2))
        !> Structural arrays on the (shared) reactive zone; the per-water transformed
        !> constants (logK_tilde/logK_star) are rebuilt here from the modified K.
        block
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
            !> Rebuild this water's per-reaction T-dependent logK FIRST (compute_logK_star consumes it
            !> via get_eq_csts_aq_chem = 10**logK), then its transformed constants logK_tilde/logK_star.
            call this%compute_logK_PHREEQC_aq_chem() !> rebuild this water's per-reaction logK at its own temperature
            call this%compute_logK_star() !> change_spec_alg uses component elimination (flag_comp=.true.): rebuild logK_star
        end block
    end subroutine

    !    !!> This subroutine sets the position of the aqueous chemistry object in the global state vector
    
    !> \brief Computes mean equilibrium reaction rate during a time step
    !!
    !! \details Calculates re_mean = Re / Delta_t for aqueous, solid,
    !! and gas chemistry.
    !!
    !! @param Delta_t time step
    !!
    subroutine compute_re_mean(this,Delta_t) !> computes mean reaction rate during a time step
        implicit none
        class(aqueous_chemistry_c) :: this !> aqueous chemistry object
        real(kind=8), intent(in) :: Delta_t !> time step
        this%rxn_rate_eq=this%rxn_amnt_eq/Delta_t
        this%solid_chemistry%rxn_rate_eq=this%solid_chemistry%rxn_amnt_eq/Delta_t
        this%gas_chemistry%rxn_rate_eq=this%gas_chemistry%rxn_amnt_eq/Delta_t
    end subroutine

    !> Compute the concentration of ideal water in the solution
    !! Mathematical Context: Water molality = 1000 g / 18.015 g/mol ≈ 55.51 mol/kg H2O
    !! For ideal dilute solutions, water activity is approximately 1, and this concentration is used for activity calculations
    subroutine compute_conc_ideal_water(this)
        !!> This function computes the concentration (in molality) of ideal water
        class(aqueous_chemistry_c) :: this                                 !< Aqueous chemistry object (polymorphic)
        this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=1d0/18d-3 !< Set water molality: 1000g/18g/mol ≈ 55.51 mol/kg H2O
    end subroutine

    !> Deep copy assignment of aqueous chemistry object (copy constructor pattern)
    !! Programming Context: This implements deep copy semantics for complex object with pointers and allocatable arrays
    !! All member variables, arrays, and pointer associations are properly copied to avoid shallow copy issues
    subroutine copy_aq_chem(this,other)
        !!> This subroutine assigns the content of an aqueous chemistry object to another one
        class(aqueous_chemistry_c) :: this                                 !< Destination aqueous chemistry object (target for assignment)
        class(aqueous_chemistry_c), intent(in) :: other                    !< Source aqueous chemistry object (source for copying)

        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        call this%copy_aq_chem_Lagr(other)                            !< Assign parent class data from local_chemistry_c base class
        !> full copy: unlike the Lagrangian parcel copy, the temperature and the quantities that go
        !> with it (pure-water density, Debye-Hückel A and B) are transferred too
        this%temp=other%temp
        this%density=other%density
        this%params_aq_sol=other%params_aq_sol
        if (allocated(other%conc_mix)) then
            this%conc_mix=other%conc_mix                        !< Copy concentration mixture constraints array (value assignment)
        end if
        if (allocated(other%conc_hat)) then
            this%conc_hat=other%conc_hat                        !< Copy the mixed concentrations entering the Newton speciation (value assignment)
        end if
        if (allocated(other%saturation_mins)) then
            this%saturation_mins=other%saturation_mins                     !< Copy mineral saturation indices Omega (last speciation of the parcel)
        end if
        !> Set pointer association to solid chemistry object (not deep copy). The kinetic reaction sets
        !> (indices_rk, indices_rk_old) are defined by the solid chemistry the water points to (the
        !> kinetic minerals of its mineral zone), so they are copied only together with that pointer;
        !> if the copy is later pointed to a different solid chemistry, set_solid_chemistry rebuilds them.
        if (associated(other%solid_chemistry)) then                        !< Check if solid chemistry pointer is associated in source
            call this%set_solid_chemistry(other%solid_chemistry)           !< Set pointer association to same solid chemistry object
            this%indices_rk=other%indices_rk                               !< Kinetic set of that solid chemistry (value assignment)
            this%indices_rk_old=other%indices_rk_old                       !< Previous-step kinetic set (aqueous kinetic + kinetic minerals of that solid chemistry)
        end if
        
        !> Set pointer association to gas chemistry object (not deep copy)
        if (associated(other%gas_chemistry)) then                          !< Check if gas chemistry pointer is associated in source
            call this%set_gas_chemistry(other%gas_chemistry)               !< Set pointer association to same gas chemistry object
        end if
        
        !> Set pointer association to aqueous phase object (not deep copy)
        
        !> Deep copy all index arrays (independent copies needed for different chemistry objects).
        !> The *_old arrays are set independently (update_attributes_dep_solid_chem, end of step), so
        !> they are unallocated on a freshly read water: each is copied on its own allocation status.
        if (allocated(other%ind_prim_species)) then                        !< Check if primary species indices are allocated in source
            this%ind_prim_species=other%ind_prim_species                   !< Deep copy primary species indices array
        end if
        if (allocated(other%ind_prim_species_old)) then
            this%ind_prim_species_old=other%ind_prim_species_old           !< Deep copy old primary species indices array
        end if
        if (allocated(other%ind_var_act_species)) then                     !< Check if variable activity species indices are allocated in source
            this%ind_var_act_species=other%ind_var_act_species             !< Deep copy variable activity species indices array
        end if
        if (allocated(other%ind_sec_species)) then                         !< Check if secondary species indices are allocated in source
            this%ind_sec_species=other%ind_sec_species                     !< Deep copy secondary species indices array
        end if
        if (allocated(other%ind_sec_species_old)) then
            this%ind_sec_species_old=other%ind_sec_species_old             !< Deep copy old secondary species indices array
        end if
        if (allocated(other%logK)) then                                   !< Per-water T-dependent equilibrium constants K_eq(T)
            this%logK=other%logK                                          !< Deep copy so per-water temperature dependence is preserved
        end if
        if (allocated(other%logK_tilde)) then                              !< Per-water transformed constants (flag_comp=.false.)
            this%logK_tilde=other%logK_tilde                               !< Deep copy so shared reactive zone is not relied upon
        end if
        if (allocated(other%logK_star)) then                               !< Per-water transformed constants (flag_comp=.true.)
            this%logK_star=other%logK_star                                 !< Deep copy so shared reactive zone is not relied upon
        end if
    end subroutine                                                         !< End of assignment subroutine

    !> Lagrangian parcel copy: moves the water state of `other` into the cell held by `this`.
    !> The temperature is NOT copied: it is a property of the cell (aquifer) the parcel moves into,
    !> not of the parcel, so an inflow water entering a warmer/colder aquifer takes the cell's
    !> temperature (its logK must then be rebuilt with compute_logK_PHREEQC_aq_chem). The same holds
    !> for the quantities that depend only on the cell temperature: the pure-water density and the
    !> Debye-Hückel parameters A and B (params_aq_sol), which are kept from the destination too.
    subroutine copy_aq_chem_Lagr(this,other)
        !!> This subroutine assigns the content of an aqueous chemistry object to another one
        class(aqueous_chemistry_c) :: this                                 !< Destination aqueous chemistry object (target for assignment)
        class(aqueous_chemistry_c), intent(in) :: other                    !< Source aqueous chemistry object (source for copying)
        real(kind=8) :: temp_cell                                          !< temperature of the destination cell, preserved across the copy
        real(kind=8) :: density_cell                                       !< pure-water density of the destination cell, preserved across the copy

        temp_cell=this%temp
        density_cell=this%density
        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        call this%copy_local_chemistry(other)                            !< Assign parent class data from local_chemistry_c base class
        this%temp=temp_cell                                                !< the parcel takes the cell's temperature
        this%density=density_cell                                          !< and the cell's pure-water density (function of its temperature)
        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        
        !> Copy scalar position array if allocated in source
        !if (allocated(other%pos)) then                                     !< Check if position array is allocated in source object
        
        !> Copy all scalar thermodynamic properties (params_aq_sol, the Debye-Hückel A and B at the
        !> cell temperature, is a property of the destination cell and is NOT copied)
        this%ionic_strength=other%ionic_strength                           !< Copy ionic strength value (mol/kg)
        this%pH=other%pH                                                   !< Copy pH value (-log10[H+])
        this%pe=other%pe                                                   !< Copy pe value (-log10[e-], electron activity)
        this%salinity=other%salinity                                       !< Copy salinity value (mass fraction or similar unit)
        this%alkalinity=other%alkalinity                                   !< Copy alkalinity value (charge balance, eq/L)
        
       
        !> Set pointer association to aqueous phase object (not deep copy)
        if (associated(other%aq_phase)) then                               !< Check if aqueous phase pointer is associated in source
            call this%set_aq_phase(other%aq_phase)                         !< Set pointer association to same aqueous phase object
        end if
        
        !> Deep copy all index arrays (independent copies needed for different chemistry objects)
        if (allocated(other%indices_aq_species)) then                      !< Check if aqueous species indices are allocated in source
            this%indices_aq_species=other%indices_aq_species               !< Deep copy aqueous species indices array
        end if
        if (allocated(other%indices_aq_phase)) then                        !< Check if aqueous phase indices are allocated in source
            this%indices_aq_phase=other%indices_aq_phase                   !< Deep copy aqueous phase indices array
        end if
        !> the parcel's previous-step state was equilibrated in the SOURCE cell, so its old solid
        !> chemistry (old reactive zone) is the source's and its old AQUEOUS kinetic rates travel with
        !> it (copied with the local chemistry above). Its old KINETIC SET (indices_rk_old) is NOT
        !> taken from the source: the old rates of the kinetic MINERALS belong to the cell the parcel
        !> is in now (solid_chemistry%rxn_rate_kin_old, see get_rk_old), so indices_rk_old keeps the
        !> destination's set (aqueous kinetic reactions + this cell's kinetic minerals, as set by
        !> update_state_variables)
        call this%set_solid_chemistry_old(other%solid_chemistry)
        !> and the previous-step index arrays are the SOURCE's current ones: the layout its
        !> previous-step concentrations (conc_old, copied above) are read with in get_conc_old_var_act_new_rz
        if (allocated(other%ind_prim_species)) this%ind_prim_species_old=other%ind_prim_species
        if (allocated(other%ind_sec_species)) this%ind_sec_species_old=other%ind_sec_species
    end subroutine    
    
        !> \brief Sets indices of aqueous variable activity species in the state variables
        ! subroutine set_ind_aq_var_act_species(this)
        !     implicit none
        !     class(aqueous_chemistry_c) :: this !> aqueous chemistry object

        !     integer(kind=4) :: i !> loop indices

        !     if (allocated(this%ind_aq_var_act_species)) then
        !         deallocate(this%ind_aq_var_act_species) !> deallocate if already allocated
        !     end if

        !     allocate(this%ind_aq_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !> allocate array for target waters

        !     do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> for each aqueous primary species
        !         this%ind_aq_var_act_species(i)=this%ind_var_act_species(i) !> set index
        !     end do
        !     do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !> for each aqueous secondary variable activity species
        !         this%ind_aq_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)=&
        !             this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i) !> set index
        !     end do
        ! end subroutine
        
        !> \brief Computes logarithmic activity coefficients for the aqueous phase
        !!
        !! \details Calculates ionic strength and then uses aqueous solution
        !! parameters to compute log activity coefficients for all aqueous species.
        !!
        subroutine compute_log_act_coeffs_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: i,k !> dissolved solid counter, its index in the aqueous phase
            real(kind=8), allocatable :: log_act_coeffs(:) !> log activity coefficients in AQUEOUS-PHASE order (compute_log_act_coeffs_aq_phase output)
            allocate(log_act_coeffs(this%aq_phase%num_species))
            call this%compute_ionic_strength()
            log_act_coeffs = 0d0 !> compute_log_act_coeffs_aq_phase sets the dissolved solids only (water is constant activity)
            call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_strength,this%params_aq_sol,log_act_coeffs)
            !> aqueous-phase order -> this water's concentration order (this%log_act_coeffs is read by
            !> concentration index everywhere: get_log_gamma1/2v, compute_activities_diss_solids). Only the
            !> dissolved solids: the water's own log gamma is set by compute_log_act_coeff_wat and kept here.
            do i=1,this%aq_phase%num_species-this%aq_phase%wat_flag
                k=this%aq_phase%ind_diss_solids(i)
                this%log_act_coeffs(this%indices_aq_species(k))=log_act_coeffs(k)
            end do
            deallocate(log_act_coeffs)
        end subroutine

        !> \brief Computes the Jacobian of logarithmic activity coefficients
        !!
        !! \details Computes the Jacobian matrix of log activity coefficients
        !! with respect to concentrations using derivatives of the activity
        !! coefficient model with respect to ionic strength.
        !!
        subroutine compute_log_Jacobian_act_coeffs_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            
            integer(kind=4) :: n_v_aq
            real(kind=8), allocatable :: d_log_gamma_d_I(:),out_prod(:,:)
            allocate(out_prod(n_v_aq,n_v_aq))
            d_log_gamma_d_I=this%compute_d_log_gamma_d_I_aq_chem() !> Compute Ã¢Ë†â€š(log ÃŽÂ³Ã¡ÂµÂ¢)/Ã¢Ë†â€šI for all variable activity species (from Debye-HÃƒÂ¼ckel or extended models)
            out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:n_v_aq)) !> Compute outer product: [Ã¢Ë†â€š(log ÃŽÂ³)/Ã¢Ë†â€šI] Ã¢Å â€” [zÃ‚Â²] for activity corrections in Jacobian
            call this%aq_phase%compute_log_Jacobian_act_coeffs_aq_phase(out_prod,this%concentrations(this%ind_var_act_species),&
                this%log_Jacobian_act_coeffs(1:n_v_aq,1:n_v_aq))
            deallocate(d_log_gamma_d_I)
        end subroutine

        !> \brief Get logarithmic activity coefficients of all secondary species.
        !>
        !> This function retrieves the logarithmic activity coefficients (log gamma) for all secondary species
        !> in the chemical system. Secondary species include aqueous secondary species, equilibrium minerals
        !> (both variable and constant activity), surface complexes (ion exchange), and equilibrium gases.
        !>
        !> The log activity coefficients are assembled in the following order:
        !> 1. Aqueous secondary species (from aqueous phase)
        !> 2. Equilibrium variable activity minerals
        !> 3. Equilibrium constant activity minerals
        !> 4. Surface complexes (cation exchange sites)
        !> 5. Equilibrium variable activity gases
        !> 6. Equilibrium constant activity gases
        !>
        !> Note: Kinetic minerals and gases are excluded from secondary species (handled separately)
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state
        !> \param[out] log_gamma2 Logarithmic activity coefficients of all secondary species (allocated internally)
        !> \return log_gamma2 Vector of log gamma values [-] (size = num_sec_species)
        function get_log_gamma2(this) result(log_gamma2)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma2(:)
            
            integer(kind=4) :: i, n2

            allocate(log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species))
            n2=0 !> counter for secondary species
            !> Aqueous secondary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                log_gamma2(n2+i)=this%log_act_coeffs(this%ind_sec_species(i))
            end do
            n2=n2+this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                log_gamma2(n2+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            n2=n2+this%solid_chemistry%reactive_zone%num_minerals_var_act
            !> Equilibrium constant activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
                log_gamma2(n2+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)
            end do
            n2=n2+this%solid_chemistry%reactive_zone%num_minerals_cst_act
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma2(n2+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            n2=n2+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma2(n2+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n2=n2+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            !> Equilibrium constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                log_gamma2(n2+i)=this%gas_chemistry%log_act_coeffs(i)
            end do
            n2=n2+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        end function

        !> \brief Get logarithmic activity coefficients of all primary species.
        !>
        !> This function retrieves the logarithmic activity coefficients (logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³Ã¢â€šÂ) for all primary
        !> species in the chemical system. Primary species (also called basis species or components)
        !> are the fundamental independent species used as a basis for expressing all other species
        !> through mass action relationships in chemical equilibrium calculations.
        !>
        !> The log activity coefficients are assembled in the following order:
        !> 1. Aqueous primary species (basis species for aqueous complexation, e.g., CaÃ‚Â²Ã¢ÂÂº, HCOÃ¢â€šÆ’Ã¢ÂÂ»)
        !> 2. Free sites (unoccupied surface sites for cation exchange, if present)
        !> 3. Kinetic variable activity minerals (minerals undergoing kinetic dissolution/precipitation)
        !> 4. Kinetic variable activity gases (gases with time-dependent partial pressures)
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state and activity coefficients
        !> \param[out] log_gamma1 Logarithmic activity coefficients of all primary species (allocated internally)
        !> \return log_gamma1 Vector of logÃ¢â€šÂÃ¢â€šâ‚¬ ÃŽÂ³Ã¢â€šÂ values [-] (size = num_prim_species)
        !>
        !> \note Primary species are independent variables in speciation calculations. Their
        !>       concentrations/activities are either specified or solved for, while secondary
        !>       species are computed from primaries via mass action law.
        !>
        !> \warning Some code at the end of this function references undefined variables (num_sp, log_gamma_vec)
        !>          and appears to be dead code that should be removed or fixed.
        !>
        !> \see get_log_gamma2, get_a1, get_c1, compute_log_act_coeffs_aq_chem
        function get_log_gamma1(this) result(log_gamma1) !> gets log_10 activity coefficients of primary species
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma1(:)

            integer(kind=4) :: i, n_p !> loop index and counter

            allocate(log_gamma1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            n_p=0 !> counter for primary species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_gamma1(i)=this%log_act_coeffs(this%ind_prim_species(i))
            end do
            n_p=n_p+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Free site
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma1(n_p+1)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
            end do
            n_p=n_p+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Kinetic variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
                log_gamma1(n_p+i)=this%solid_chemistry%log_act_coeffs(i)
            end do
            n_p=n_p+this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
            !> Kinetic variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
                log_gamma1(n_p+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+i)
            end do
            n_p=n_p+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
        end function get_log_gamma1

        function get_log_gamma_var_act(this) result(log_gamma_var_act)
        !!> This function returns the logarithmic activity coefficients of all variable activity species
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma_var_act(:)

            integer(kind=4) :: i, num_sp

            allocate(log_gamma_var_act(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            num_sp=0 !> counter for variable activity species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_gamma_var_act(i)=this%log_act_coeffs(this%ind_var_act_species(i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma_var_act(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
                log_gamma_var_act(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)=&
                    this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Kinetic variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
                log_gamma_var_act(num_sp+i)=this%solid_chemistry%log_act_coeffs(i)
            end do
            num_sp=num_sp+this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
            !> Kinetic variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
                log_gamma_var_act(num_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
            !> Aqueous secondary variable activity
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                log_gamma_var_act(num_sp+i)=this%log_act_coeffs(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                log_gamma_var_act(num_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%num_minerals_var_act+&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats !> chapuza
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma_var_act(num_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            num_sp=num_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            !!> Ideal water
            !!> Kinetic constant activity minerals
            !!> Equilibrium constant activity minerals
            !!> Kinetic constant activity gases
            !!> Equilibrium constant activity gases
        end function get_log_gamma_var_act

        !> \brief Computes logarithmic activity coefficients assuming ideal solution
        !!
        !! \details Sets all log activity coefficients to zero (ideal solution
        !! assumption: gamma = 1 for all variable activity species).
        !!
        !! @return log_gamma_var_act vector of log activity coefficients (all zeros for ideal)
        !!
        function compute_log_act_coeffs_var_act_ideal(this) result(log_gamma_var_act)
            implicit none
            ! Local variables
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            real(kind=8), allocatable :: log_gamma_var_act(:)
            
            integer(kind=4) :: i, n_v
            
            allocate(log_gamma_var_act(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            n_v=0 !> counter for variable activity species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_gamma_var_act(i)=0d0 !> ideal solution: log ÃŽÂ³ = 0
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma_var_act(n_v+1)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
                log_gamma_var_act(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)=&
                    this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Kinetic variable activity minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
                log_gamma_var_act(n_v+i)=this%solid_chemistry%log_act_coeffs(i)
            end do
            n_v=n_v+this%solid_chemistry%mineral_zone%num_minerals_kin_var_act
            !> Kinetic variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
                log_gamma_var_act(n_v+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+i)
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act
            !> Aqueous secondary variable activity
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                log_gamma_var_act(n_v+i)=0d0 !> ideal solution: log ÃŽÂ³ = 0
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                log_gamma_var_act(n_v+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%num_minerals_var_act+&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats !> chapuza
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma_var_act(n_v+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n_v=n_v+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
        end function compute_log_act_coeffs_var_act_ideal

         !> \brief Computes logarithmic activity coefficients assuming ideal solution
        !!
        !! \details Sets all log activity coefficients to zero (ideal solution
        !! assumption: gamma = 1 for all variable activity species).
        !!
        !! @return log_gamma_ideal vector of log activity coefficients (all zeros for ideal)
        !!
        function compute_log_act_coeffs_ideal(this) result(log_gamma_ideal)
            implicit none
            ! Local variables
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            real(kind=8), allocatable :: log_gamma_ideal(:)
            
            integer(kind=4) :: i, n_sp
            
            allocate(log_gamma_ideal(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
            n_sp=0 !> counter for total species
            !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_gamma_ideal(i)=0d0 !> ideal solution: log gamma = 0
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            !> Free site
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma_ideal(n_sp+1)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1) !> free site (chapuza)
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl-&
                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Aqueous secondary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                log_gamma_ideal(n_sp+i)=0d0 !> ideal solution: log gamma = 0
            end do            
            n_sp=n_sp+this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            !> Equilibrium variable activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                log_gamma_ideal(n_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%num_minerals_var_act
            !> Equilibrium constant activity minerals
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
                log_gamma_ideal(n_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+i)
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%num_minerals_cst_act
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                log_gamma_ideal(n_sp+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i) !> surface complexes
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Equilibrium variable activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                log_gamma_ideal(n_sp+i)=this%gas_chemistry%log_act_coeffs(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
            !> Equilibrium constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                log_gamma_ideal(n_sp+i)=this%gas_chemistry%log_act_coeffs(i)
            end do
            n_sp=n_sp+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        end function compute_log_act_coeffs_ideal
        
!> \brief Positivity-preserving AND component-conservative reactive update of the mixed
!!        concentrations (conservative reaction-extent limiter).
!!
!! \details Replaces the explicit reactive update c_hat = c_mix + Delta_t*r_tilde, which can
!! produce negative concentrations when the (explicit) reaction contribution r_tilde acts as a
!! strong sink, WITHOUT breaking the component mass balance the WMA relies on. Downstream the
!! solver uses u_hat = U*c_hat as the conserved component total (see reactive_mixing_iter_*),
!! so the positivity fix must preserve U*c_hat, not merely keep c_hat >= 0.
!!
!! \par Why not (classical or modified) Patankar
!! The classical per-species Patankar update (c_hat_i = c_mix_i^2/(c_mix_i - Delta_t*r_tilde_i)
!! on sinks) is unconditionally positive but NOT conservative: it damps only the negative
!! entries of r_tilde, so U*c_hat drifts from U*(c_mix + Delta_t*r_tilde). The textbook
!! "modified" Patankar (Burchard, Deleersnijder & Meister 2003) restores conservation of the
!! PLAIN SUM sum(c) only; it does NOT conserve the component totals u = U*c that are the actual
!! transported invariants here (and sum(c) is not even reaction-invariant for non-unit
!! stoichiometry). Moreover r_tilde is an already-MIXED net vector (sum_j lambda_j S_j^T r_j over
!! several waters/zones), so it carries no single clean per-reaction production-destruction
!! structure to build a pairwise coupled system from.
!!
!! \par Method (reaction-extent limiter)
!! Scale the WHOLE net reaction vector by a single factor omega in [0,1]:
!!   c_hat = c_mix + Delta_t*omega*r_tilde,
!!   omega = min(1, min_{i: r_tilde_i<0} c_mix_i / (Delta_t*(-r_tilde_i))).
!! Because omega multiplies r_tilde as a whole, the update stays parallel to r_tilde, hence
!!   U*c_hat = U*c_mix + Delta_t*omega*(U*r_tilde).
!! Equilibrium reactions satisfy U*Se^T = 0, so U*r_tilde has NO equilibrium part: the
!! equilibrium redistribution (the usual source of negative species) is conserved EXACTLY
!! (U*c_hat = U*c_mix) for any omega, while the kinetic component change stays a consistent,
!! omega-sub-stepped valid reaction (it remains in range(U*Sk^T), so no spurious components are
!! created). Positivity is unconditional for c_mix >= 0 and any Delta_t. The reaction mixing
!! ratios (and hence the discrete dispersion) are left untouched, unlike modify_mix_ratios_reacts.
!!
!! \note omega is global: a single near-empty species being consumed throttles the whole reaction
!! step. For stiff multi-kinetic problems a per-reaction (reaction-resolved) limiter is more
!! accurate, but it needs un-mixed stoichiometry/rate data unavailable at this post-mixing stage.
!!
!! \par References
!! - Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow. Hemisphere, New York.
!! - Burchard, H., Deleersnijder, E. & Meister, A. (2003). A high-order conservative
!!   Patankar-type discretisation for stiff systems of production-destruction equations.
!!   Applied Numerical Mathematics 47(1), 1-30.
!!
!! @param c_mix mixed concentrations after transport (assumed >= 0)
!! @param Delta_t time step
!! @param r_tilde net (already mixed) reaction rate contributions (same ordering as c_mix)
!! @param c_hat positivity-preserving, component-conservative updated concentrations (already allocated)
subroutine compute_c_hat_patankar(c_mix,Delta_t,r_tilde,c_hat)
    implicit none
    real(kind=8), intent(in) :: c_mix(:)   !> mixed concentrations after transport
    real(kind=8), intent(in) :: Delta_t    !> time step
    real(kind=8), intent(in) :: r_tilde(:) !> net reaction rate contributions after mixing
    real(kind=8), intent(out) :: c_hat(:)  !> updated concentrations (must be already allocated)
    integer(kind=4) :: i           !> species index
    real(kind=8) :: omega          !> global reaction-extent factor in [0,1]
    real(kind=8) :: omega_i        !> per-species positivity bound
    !> Largest omega in [0,1] keeping every species non-negative under a whole-vector scaling.
    omega=1d0
    do i=1,size(c_mix)
        if (r_tilde(i)<0d0) then
            !> denom = Delta_t*(-r_tilde(i)) > 0 for Delta_t>0; omega_i = c_mix(i)/denom >= 0.
            omega_i=c_mix(i)/(Delta_t*(-r_tilde(i)))
            if (omega_i<omega) omega=omega_i
        end if
    end do
    if (omega<0d0) omega=0d0
    !> Whole-vector scaling: preserves U*c_hat exactly for equilibrium redistribution
    !> (U*r_tilde has no equilibrium part) and sub-steps the kinetic part consistently.
    c_hat=c_mix+Delta_t*omega*r_tilde
end subroutine

!> Modify mixing ratios for kinetic reaction rates to prevent negative concentrations
!! Mathematical Context: Adaptive time stepping for reactive transport to maintain positivity
!! Algorithm: If c + ÃŽâ€tÃƒâ€”r < 0, reduce mixing ratio ÃŽÂ» to ensure c + ÃŽâ€tÃƒâ€”ÃŽÂ»Ãƒâ€”r Ã¢â€°Â¥ 0
!! Physical Context: Prevents unphysical negative concentrations during reactive mixing
subroutine modify_mix_ratios_reacts(this,mix_ratios_old,c_mix,Delta_t,r_tilde,&
    mix_ratios_new,num_lump)
    !> This subroutine modifies the mixing ratios of the reaction amounts for a target
    !> AQUI DEBERIAS GUARDAR LOS NUEVOS LAMBDAS                           !< TODO: Should save new lambdas
        class(aqueous_chemistry_c), intent(in) :: this                       !< Solid chemistry object (input only)
        real(kind=8), intent(in) :: mix_ratios_old(:)                         !< Initial mixing ratio of reaction rates
        real(kind=8), intent(in) :: c_mix(:)                             !< Concentration vector after transport
        real(kind=8), intent(in) :: Delta_t                                !< Time step size
        !real(kind=8), intent(in) :: theta                                   !< Time weighting factor
        real(kind=8), intent(inout) :: r_tilde(:)                         !< reaction rate contributions
        real(kind=8), intent(inout) :: mix_ratios_new(:)                      !< New mixing ratio for kinetic reactions (output)
        integer(kind=4), intent(out) :: num_lump                           !< Number of lumping operations performed

        integer(kind=4) :: i                                               !< Loop index for species
        integer(kind=4) :: n_v                                            !< Number of variable activity species
        real(kind=8), parameter :: alpha=1.1                              !< Mixing ratio modification factor (5% increase)
        real(kind=8) :: beta                              !< Mixing ratio modification factor (5% increase)
        real(kind=8), allocatable :: rk_old(:)                            !< Reaction amount contributions after mixing [mol/L]
        real(kind=8), allocatable :: r_tilde_old(:)                            !< Reaction amount contributions after mixing [mol/L]
        real(kind=8), allocatable :: Sk_nc(:,:)                           !< Kinetic reaction stoichiometric matrix for non-conserved species
        real(kind=8) :: lumped_lambda                !< New mixing ratios for kinetic reactions

        rk_old=this%get_rk_old_new_set()                           !< Previous-step kinetic rates in the current kinetic set (rows of Sk_nc)
        Sk_nc=this%get_Sk_nc()
        n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species          !< Get number of variable activity species
        mix_ratios_new=mix_ratios_old                                 !< Initialize new mixing ratio (temporary - "chapuza")
        lumped_lambda=sum(mix_ratios_old)                          !< Initialize lumped mixing ratio
        r_tilde_old=r_tilde                                                     !< Commented lumping flag initialization
        num_lump=0                                                          !< Initialize number of lumping operations
        i=1                                                                 !< Initialize species index
        !alpha=1d0
        beta=1d0
        !> Check for negative concentrations after mixing and adjust if necessary
        do                                                                  !< Infinite loop with exit condition
            if (c_mix(i)+Delta_t*r_tilde(i)<=0d0) then                           !< Check if total concentration would become negative
                beta=(lumped_lambda-alpha*mix_ratios_new(1))/(&
                    lumped_lambda-mix_ratios_new(1))                                             !< Commented lumping flag setting
                num_lump=num_lump+1                                        !< Increment lumping counter
                mix_ratios_new(1)=mix_ratios_new(1)*alpha                    !< Increase mixing ratio by factor ÃŽÂ±
                mix_ratios_new(2:)=mix_ratios_new(2:)*beta
                if (mix_ratios_new(1)>lumped_lambda) then                     !< Check for excessive mixing ratio
                    !> Limit mixing ratio to prevent instability (hack)
                    mix_ratios_new(1)=lumped_lambda                           !< Set maximum stable mixing ratio
                    mix_ratios_new(2:)=0d0                           !< Set maximum stable mixing ratio
                    !beta=0d0                           !< Set maximum stable mixing ratio
                    r_tilde=0d0 !< apply final mixing ratio adjustment
                    exit
                else
                    r_tilde=r_tilde*beta !< apply final mixing ratio adjustment
                end if
            else if (i<n_v) then
                i=i+1
            else
                ! mix_ratios_new(1)=mix_ratios_init(1)*alpha                    !< Increase mixing ratio by factor ÃŽÂ±
                !mix_ratios_new(2:)=mix_ratios_new(2:)*beta                    !< Increase mixing ratio by factor ÃŽÂ±
                exit
            end if
        end do
        if (abs(sum(mix_ratios_new)-lumped_lambda) > this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
            print *, sum(mix_ratios_new), lumped_lambda
            error stop "Error in mixing ratios conservation in subroutine modify_mix_ratios_reacts"
        end if
end subroutine

subroutine compute_r_tilde_lump(this,c_mix,Delta_t,theta,lambda_r,r_tilde)
    !> This subroutine computes the old reaction rate contributions in the lumped scheme
        class(aqueous_chemistry_c), intent(in) :: this                       !< Solid chemistry object (input only)
        real(kind=8), intent(in) :: c_mix(:)                             !< Concentration vector after transport
        real(kind=8), intent(in) :: Delta_t                                !< Time step size
        real(kind=8), intent(in) :: theta                                   !< Time weighting factor
        !real(kind=8), intent(in) :: rk_old(:)                             !< Old reaction amount contributions [mol/L]
        !real(kind=8), intent(in) :: Sk_nc(:,:)                           !< Kinetic reaction stoichiometric matrix for non-conserved species
        real(kind=8), intent(inout) :: lambda_r                      !< New mixing ratio for kinetic reactions
        !class(aqueous_chemistry_c), intent(in) :: upstream_water                       !< Solid chemistry object (input only)
        real(kind=8), intent(out) :: r_tilde(:)                         !< Modified reaction amount contributions [mol/L]

        integer(kind=4) :: i                                               !< Loop index for species
        real(kind=8), parameter :: alpha=1.05                              !< Mixing ratio modification factor (5% decrease)
        real(kind=8), allocatable :: rk_old(:)                            !< Reaction rates after previous time step
        real(kind=8), allocatable :: Sk_nc(:,:)                           !< Kinetic stoichiometric matrix for non-conserved species
        !> Compute modified reaction amount contributions after mixing
        rk_old=this%get_rk_old_new_set()                           !< Previous-step kinetic rates in the current kinetic set (rows of Sk_nc)
        Sk_nc=this%get_Sk_nc()                        !< Get kinetic stoichiometric matrix for non-conserved species
        r_tilde=(1d0-theta)*lambda_r*matmul(transpose(Sk_nc),rk_old)
        i=1
        do
            if (c_mix(i)+Delta_t*r_tilde(i)<0d0) then
                lambda_r=lambda_r*alpha !< reduce mixing ratio
                if (lambda_r>1d0) then
                    lambda_r=1d0
                    r_tilde=(1d0-theta)*matmul(transpose(Sk_nc),rk_old)
                    exit
                end if
                r_tilde=r_tilde*alpha !< we should not reach this point
                !error stop "Negative concentration after mixing in subroutine compute_r_tilde_lump"
            else if (i<size(r_tilde)) then
                i=i+1
            else
                exit
            end if
        end do
        deallocate(rk_old)
        deallocate(Sk_nc)
end subroutine

subroutine update_solid_chemistry_old(this)
    !> This subroutine updates the solid chemistry old state to the current state
        class(aqueous_chemistry_c), intent(inout) :: this                       !< Solid chemistry object (input/output)
        this%solid_chemistry_old=>this%solid_chemistry !< Update old solid chemistry state
        !call this%solid_chemistry%update_old_state()
end subroutine

    subroutine update_indices_rk_old(this)
    !> This subroutine updates the indices for retrieving old reaction rates
        class(aqueous_chemistry_c), intent(inout) :: this                       !< Solid chemistry object (input/output)
        ! integer(kind=4) :: i                                               !< Loop index
        ! integer(kind=4) :: n_k                                            !< Counter for kinetic reactions
        this%indices_rk_old=this%indices_rk !< Update indices for old reaction rates
    end subroutine update_indices_rk_old

    function get_re_mean(this) result(re_mean)
        !> This function gets the mean equilibrium reaction rates during a time step
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: re_mean(:)

            integer(kind=4) :: i,n_eq !> loop index and counter

            allocate(re_mean(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            
            n_eq=0 !> counter for equilibrium reactions
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_cst_act
                re_mean(n_eq+i)=this%solid_chemistry%rxn_rate_eq(i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%num_minerals_cst_act
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                re_mean(n_eq+i)=this%gas_chemistry%rxn_rate_eq(i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
                re_mean(n_eq+i)=this%rxn_rate_eq(n_eq+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                re_mean(n_eq+i)=this%solid_chemistry%rxn_rate_eq(&
                    this%solid_chemistry%reactive_zone%num_minerals_cst_act+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%num_minerals_var_act
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                re_mean(n_eq+i)=this%solid_chemistry%rxn_rate_eq(&
                    this%solid_chemistry%reactive_zone%num_minerals+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                re_mean(n_eq+i)=this%gas_chemistry%rxn_rate_eq(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
        end function get_re_mean

        subroutine update_state_variables(this)
        !> This subroutine updates the state variables of the aqueous & solid chemistry objects
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            call this%update_conc_old() !< Update the aqueous chemistry state variables
            call this%solid_chemistry%update_conc_old() !< Update the solid chemistry state variables
            !call this%gas_chemistry%update_conc_old() !< Update the gas chemistry state variables
            call this%update_rxn_rates_and_amounts() !< Update the aqueous reaction rates
            call this%solid_chemistry%update_rxn_rates_and_amounts() !< Update the solid reaction rates
            
        end subroutine update_state_variables

        subroutine update_attributes_dep_solid_chem(this)
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            !> This subroutine updates the attributes that depend on the solid chemistry associated with the aqueous chemistry object
            this%indices_rk_old=this%indices_rk !< Update indices for old reaction rates
            this%solid_chemistry_old=>this%solid_chemistry !< Update old solid chemistry state
            this%ind_prim_species_old=this%ind_prim_species !< Update old primary species indices
            this%ind_sec_species_old=this%ind_sec_species !< Update old secondary species indices
        end subroutine update_attributes_dep_solid_chem

        !> \brief Shifts the time levels of the water at the start of a step: current -> old (and
        !!        old -> old_old), as update_old_attributes did before the split.
        !! \details 1. state (update_state_variables): aqueous and solid concentrations, reaction rates
        !!             and amounts; 2. the attributes that depend on the solid chemistry
        !!             (update_attributes_dep_solid_chem): previous-step kinetic set, solid chemistry and
        !!             index arrays; 3. the attributes that depend on the (possibly new) reactive zone
        !!             (set_attributes_dep_rz): index arrays, logK, logK_star. Without step 1 the mixers
        !!             read a stale conc_old and re-mix the initial state every step.
        subroutine update_attributes(this) !< Update all attributes of the aqueous chemistry object
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            call this%update_state_variables() !< Shift the state: concentrations, reaction rates and amounts (current -> old)
            call this%update_attributes_dep_solid_chem() !< Update attributes dependent on solid chemistry
            !call this%set_attributes_dep_rz() !< Set attributes dependent on the reactive zone
        end subroutine update_attributes !< Update all attributes of the aqueous chemistry object
        

        subroutine set_attributes_dep_rz(this)
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            !> This subroutine sets the attributes that depend on the reactive zone associated with the aqueous chemistry object
            call this%set_ind_species_react_zone()
            call this%compute_logK_PHREEQC_aq_chem()
            call this%compute_logK_star()
        end subroutine set_attributes_dep_rz

        subroutine set_attributes_dep_mz(this)
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            !> This subroutine sets the attributes that depend on the mineral zone associated with the aqueous chemistry object
            call this%compute_saturation_mins()
            call this%set_indices_rk()
        end subroutine set_attributes_dep_mz

        !> \brief Set the mean equilibrium reaction rates during a time step
        !> \param this Local chemistry object
        !> \param Delta_t Time step size
        subroutine set_eq_rxn_rates(this,Delta_t)
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: Delta_t
            this%rxn_rate_eq = this%rxn_amnt_eq / Delta_t !> current-step mean equilibrium rate (the rxn_rate_eq_old/rxn_rate_eq_old_old history is shifted end-of-step in update_rxn_rates)
            this%solid_chemistry%rxn_rate_eq = this%solid_chemistry%rxn_amnt_eq / Delta_t !> current-step mean equilibrium rate for solid chemistry
            if (associated(this%gas_chemistry)) then
                this%gas_chemistry%rxn_rate_eq = this%gas_chemistry%rxn_amnt_eq / Delta_t !> current-step mean equilibrium rate for gas chemistry
            end if
        end subroutine set_eq_rxn_rates

        subroutine compute_rxn_rate_kin_avg(this,theta_r)
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: theta_r
            this%rxn_rate_kin_avg = theta_r * this%rxn_rate_kin + (1d0-theta_r) * this%rxn_rate_kin_old
            this%solid_chemistry%rxn_rate_kin_avg = theta_r * this%solid_chemistry%rxn_rate_kin + (1d0-theta_r) * this%solid_chemistry%rxn_rate_kin_old
        end subroutine compute_rxn_rate_kin_avg

        subroutine set_tar_wat_idx(this,tar_wat_idx)
            class(zone_water_c), intent(inout) :: this
            integer, intent(in) :: tar_wat_idx
            if (tar_wat_idx <= 0) then
                error stop "Target water index must be positive."
            else
                this%tar_wat_idx = tar_wat_idx
            end if
        end subroutine set_tar_wat_idx

        subroutine set_zone_id(this,zone_id)
            class(zone_water_c), intent(inout) :: this
            integer, intent(in) :: zone_id
            if (zone_id <= 0) then
                error stop "Zone ID must be positive."
            else
                this%zone_id = zone_id
            end if
        end subroutine set_zone_id

        subroutine set_lambda_c(this,lambda_c)
            class(zone_water_c), intent(inout) :: this
            real(kind=8), intent(in) :: lambda_c
            if (lambda_c < 0d0) then
                error stop "Lambda_c must be non-negative."
            end if
            this%lambda_c = lambda_c
        end subroutine set_lambda_c

        subroutine set_aq_chem(this,aq_chem)
            class(zone_water_c), intent(inout) :: this
            class(aqueous_chemistry_c), intent(in) :: aq_chem
            if (.not. associated(aq_chem%solid_chemistry)) then
                error stop "Aqueous chemistry must have an associated solid chemistry."
            else if (.not. associated(aq_chem%aq_phase)) then
                error stop "Aqueous chemistry must have an associated aqueous phase."
            end if
            this%aq_chem = aq_chem
        end subroutine set_aq_chem

        !> Copies the AQUEOUS state of other into this (scalars, concentration/activity/reaction arrays,
        !> species index arrays) and points this to other's aqueous phase; the solid/gas chemistry
        !> pointers of this are left untouched (a zone water keeps its own scratch solid chemistry).
        !> Array components are copied with copy_alloc: unallocated in other => deallocated in this.
        subroutine copy_aqueous_attributes(this,other)
            class(aqueous_chemistry_c), intent(inout) :: this
            class(aqueous_chemistry_c), intent(in) :: other
            this%params_aq_sol = other%params_aq_sol
            this%pH = other%pH
            this%pe = other%pe
            this%temp = other%temp
            this%density = other%density
            this%volume = other%volume
            this%pressure = other%pressure
            this%ionic_strength = other%ionic_strength
            this%salinity = other%salinity
            this%alkalinity = other%alkalinity
            call copy_alloc(this%log_act_coeffs, other%log_act_coeffs)
            call copy_alloc(this%log_Jacobian_act_coeffs, other%log_Jacobian_act_coeffs)
            call copy_alloc(this%concentrations, other%concentrations)
            call copy_alloc(this%conc_old, other%conc_old)
            call copy_alloc(this%conc_old_old, other%conc_old_old)
            call copy_alloc(this%conc_mix, other%conc_mix)
            call copy_alloc(this%conc_hat, other%conc_hat)
            call copy_alloc(this%activities, other%activities)
            call copy_alloc(this%rxn_rate_kin, other%rxn_rate_kin)
            call copy_alloc(this%rxn_amnt_kin, other%rxn_amnt_kin)
            call copy_alloc(this%cum_rxn_amnt_kin, other%cum_rxn_amnt_kin)
            call copy_alloc(this%cum_rxn_amnt_eq, other%cum_rxn_amnt_eq)
            call copy_alloc(this%rxn_rate_kin_avg, other%rxn_rate_kin_avg)
            call copy_alloc(this%rxn_rate_kin_old, other%rxn_rate_kin_old)
            call copy_alloc(this%rk_old_old, other%rk_old_old)
            call copy_alloc(this%rk_old_old_old, other%rk_old_old_old)
            call copy_alloc(this%rxn_amnt_eq, other%rxn_amnt_eq)
            call copy_alloc(this%rxn_rate_eq, other%rxn_rate_eq)
            call copy_alloc(this%rxn_rate_eq_old, other%rxn_rate_eq_old)
            call copy_alloc(this%rxn_rate_eq_old_old, other%rxn_rate_eq_old_old)
            call copy_alloc(this%rxn_rate_eq_old_old_old, other%rxn_rate_eq_old_old_old)
            call copy_alloc(this%contrib_rxn_amnt_eq, other%contrib_rxn_amnt_eq)
            call copy_alloc(this%contrib_rxn_amnt_kin, other%contrib_rxn_amnt_kin)
            call copy_alloc(this%var_act_species_indices, other%var_act_species_indices)
            call copy_alloc(this%cst_act_species_indices, other%cst_act_species_indices)
            call copy_alloc(this%indices_aq_species, other%indices_aq_species)
            call copy_alloc(this%indices_aq_phase, other%indices_aq_phase)
            if (associated(other%aq_phase)) then
                call this%set_aq_phase(other%aq_phase)
            else
                error stop "Other aqueous chemistry object must have an associated aqueous phase."
            end if
        end subroutine copy_aqueous_attributes

        subroutine sum_eq_rxn_rates(this,zone_waters)
            class(aqueous_chemistry_c), intent(inout) :: this
            !real(dp), intent(in) :: lambda_r
            class(zone_water_c), intent(in) :: zone_waters(:)
            integer :: i
            do i=2,size(zone_waters)-1
                if ((this%solid_chemistry%is_interface .eqv. .true.) .and. this%solid_chemistry%reactive_zone%id==&
                    zone_waters(i)%zone_id) then
                        this%solid_chemistry%rxn_rate_eq = this%solid_chemistry%rxn_rate_eq + &
                            zone_waters(i-1)%aq_chem%solid_chemistry%rxn_rate_eq + &
                            zone_waters(i+1)%aq_chem%solid_chemistry%rxn_rate_eq
                        this%solid_chemistry%rxn_amnt_eq = this%solid_chemistry%rxn_amnt_eq + &
                            zone_waters(i-1)%aq_chem%solid_chemistry%rxn_amnt_eq + &
                            zone_waters(i+1)%aq_chem%solid_chemistry%rxn_amnt_eq
                        this%solid_chemistry%cum_rxn_amnt_eq = this%solid_chemistry%cum_rxn_amnt_eq + &
                            zone_waters(i-1)%aq_chem%solid_chemistry%cum_rxn_amnt_eq + &
                            zone_waters(i+1)%aq_chem%solid_chemistry%cum_rxn_amnt_eq
                end if
            end do
            !call this%compute_contrib_rxn_amnt_eq(lambda_r)
        end subroutine sum_eq_rxn_rates
end module                                                                  !< End of aqueous_chemistry_m module
