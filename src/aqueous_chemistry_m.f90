!> \file aqueous_chemistry_m.f90
!> \brief Defines the aqueous_chemistry_c type and its associated procedures for aqueous phase chemical modeling.
!>
!> This module provides the core data structures and algorithms for simulating aqueous chemistry, including species concentrations, activities, activity coefficients, and coupling with solid and gas phases. It supports equilibrium and kinetic reactions, mixing iterations, and advanced speciation algorithms for reactive transport.
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2024
module aqueous_chemistry_m
    use aq_phase_m, only: aq_phase_c                             !< Import aqueous phase class for species and reactions
    use local_chemistry_m, only: local_chemistry_c                         !< Import base local chemistry class for inheritance
    use solid_chemistry_m, only: solid_chemistry_c !< Import solid phase chemistry classes for mineral interactions
    use gas_chemistry_m, only: gas_chemistry_c                          !< Import gas chemistry for gas-water equilibrium
    use arrays_m, only: int_array_c, diag_matrix_c, id_matrix, norm_mat_inf
    use metodos_sist_lin_m, only: LU_lin_syst
    use vectors_m, only: inf_norm_vec_real, outer_prod_vec
    use params_aq_sol_m, only: params_aq_sol_s                             !< Import aqueous solution parameter structure
    use array_ops_m, only: get_ind_int, is_int_in_1D_array                 !< Import array utility functions for index operations
    implicit none                                                           !< Require explicit variable declaration (Fortran safety)
    save   !< Preserve module variables between procedure calls      
    private !< Make all entities private by default
    public ::  initialise_iterative_method                     !< Expose iterative method initialisation subroutine
    public :: mixing_iter_comp, reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt4, mixing_iter_comp_ideal, compute_r_tilde_impl_opt1, reactive_mixing_iter_EI_eq_anal_ideal, &
        mixing_iter_species, reactive_mixing_iter_EE_kin, reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_c_mix, compute_c_mix_global, compute_r_tilde_impl_opt2, compute_r_tilde_impl_opt3, compute_r_tilde_impl_opt1_bis, &
        compute_r_tilde_impl_opt4_bis !< Expose mixing subroutines for external use
    !> \brief Aqueous chemistry class for modeling chemical reactions in the aqueous phase.
    !>
    !> This type extends local_chemistry_c and encapsulates all properties and methods required for aqueous chemical speciation, activity calculations, and coupling with solid and gas chemistry. It manages concentrations, activities, activity coefficients, and provides interfaces for mixing, reaction iteration, and property calculation.
    type, public, extends(local_chemistry_c) :: aqueous_chemistry_c
        type(params_aq_sol_s) :: params_aq_sol                             !< Aqueous solution parameters: temperature, pressure, density, etc.
        real(kind=8) :: ionic_strength                                     !< Ionic strength (I) = 0.5 * Σ(ci * zi²) for activity coefficient calculations
        real(kind=8) :: pH                                                 !< pH = -log10[H+] - measure of hydrogen ion activity (acidity)
        real(kind=8) :: pe                                                 !< pe = -log10[e-] - measure of electron activity (redox potential)
        real(kind=8) :: salinity                                           !< Total dissolved solids concentration (salinity measure)
        real(kind=8) :: alkalinity                                         !< Alkalinity = m_OH- + m_HCO3- + 2*m_CO3-2 (acid neutralizing capacity)
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
        integer(kind=4), allocatable :: ind_aq_var_act_species(:)           !< Indices of aqueous variable activity species in the global state variables
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
        procedure :: set_ind_aq_var_act_species                     !< Set indices of aqueous variable activity species in global state variables
        procedure :: set_spec_alg_aq_chem                           !< Set speciation algebra object for equilibrium calculations
        procedure :: set_ind_species                                !< Set species indices for chemical system navigation
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
        procedure :: set_Rk_est                                     !< Set estimated kinetic reaction rates for predictor-corrector methods
        procedure :: set_rk_mean                                    !< Set mean kinetic reaction rates for time stepping
        procedure :: set_rk_old                                  !< Set kinetic reaction rates from previous time step (for time averaging and predictor-corrector)
    !> Allocate procedures - Memory allocation methods for dynamic arrays
        procedure :: allocate_reaction_rates                        !< Allocate memory for kinetic reaction rate arrays
        procedure :: allocate_conc_aq_species                       !< Allocate memory for aqueous species concentration arrays
        procedure :: allocate_activities_aq_species                 !< Allocate memory for aqueous species activity arrays
        procedure :: allocate_log_act_coeffs_aq_chem                !< Allocate memory for logarithmic activity coefficient arrays
        procedure :: allocate_indices_aq_species_aq_chem            !< Allocate memory for aqueous species index arrays
        procedure :: allocate_indices_aq_phase_aq_chem              !< Allocate memory for aqueous phase index arrays
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
        procedure :: compute_salinity                               !< Compute total dissolved solids (salinity)
        procedure :: compute_alkalinity                             !< Compute alkalinity (acid neutralizing capacity)
        procedure :: compute_pH                                     !< Compute pH from hydrogen ion activity
        procedure :: compute_conc_comp                              !< Compute component concentrations for mass balance
        procedure :: compute_conc_comp_cst_act                      !< Compute component concentrations with constant activity species
        procedure :: compute_re_mean                                    !< Compute mean reaction extent for time averaging
        procedure :: compute_r_tilde_lump                                !< Compute reaction extent with spatial lumping
        procedure :: compute_Re                                     !< Compute reaction extent (progress of reactions)
        procedure :: compute_Re_kin                                  !< Compute reaction extent and kinetic reaction rates
        procedure :: compute_rk_new                                 !< Compute new kinetic reaction rates for current time step
        procedure :: compute_Jacobian_rk_anal                       !< Compute analytical Jacobian of kinetic reaction rates
        procedure :: compute_rk_Jac_rk_anal                         !< Compute kinetic rates and their analytical Jacobian
        procedure :: compute_rk_Jac_rk_incr_coeff                   !< Compute kinetic rates and Jacobian using incremental coefficients
        procedure :: compute_log_K_aq_chem                          !< Compute logarithmic equilibrium constants for aqueous reactions
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
        procedure :: compute_log_Jacobian_act_coeffs_aq_chem        !< Compute logarithmic activity coefficient of water
        procedure :: compute_saturation_kin_min                     !< Compute saturation indices for kinetic minerals
        procedure :: compute_Rk                                     !< Compute kinetic reaction rates matrix
    !> Get procedures - Accessor methods for retrieving chemical data
        procedure :: get_rk_new                                     !< Get new kinetic reaction rates (current time step)
        procedure :: get_rk_old                                     !< Get old kinetic reaction rates (previous time step)
        procedure :: get_Rk                                         !< Get kinetic reaction amounts
        procedure :: get_Sk_nc                                      !< Get kinetic stoichiometric matrix for variable activity species
        procedure :: get_S_nc                                      !< Get kinetic stoichiometric matrix for variable activity species
        procedure :: get_S_nc_old
        procedure :: get_r_old
        procedure :: get_r_new                                      !< Get kinetic stoichiometric matrix for variable activity species
        procedure :: get_Sk_nc_old                                  !< Get old kinetic stoichiometric matrix for variable activity species
        procedure :: get_c1                                         !< Get primary species concentrations vector
        procedure :: get_c1_old_old                                 !< Get primary aqueous species concentrations
        procedure :: get_a1                                         !< Get primary species activities vector
        procedure :: get_conc_nc_aq_perm                            !< Get non-constant activity aqueous species concentrations (permeable)
        procedure :: get_log_gamma_var_act                          !< Get logarithmic activity coefficients for variable activity species
        procedure :: get_c2                                         !< Get secondary species concentrations vector
        procedure :: get_c2v                                        !< Get secondary non-constant activity species concentrations
        procedure :: get_conc_nc                                    !< Get all non-constant activity species concentrations
        procedure :: get_conc_nc_old                                !< Get previous time step non-constant activity concentrations
        procedure :: get_conc                                       !< Get all species concentrations (primary + secondary)
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
        procedure :: compute_dc2_dc1                                !< General Jacobian matrix computation
        procedure :: compute_c_from_u_bin_syst                      !< Concentration from binary system molality
    !> Mixing procedures for water-rock interaction and transport coupling
        procedure :: mixing_iter_comp                               !< Iterative component mixing with activity corrections
        procedure :: mixing_iter_comp_ideal                         !< Component mixing assuming ideal solutions
        procedure :: mixing_iter_species                            !< Species-level iterative mixing procedure
        procedure :: reactive_mixing_iter_EE_eq_kin                    !< Water mixing with explicit equilibrium & kinetics
        procedure :: reactive_mixing_iter_EE_eq_kin_ideal              !< Water mixing EE method (ideal solution)
        procedure :: reactive_mixing_iter_EE_kin               !< EE kinetics-only water mixing (ideal)
        procedure :: reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2    !< EI analytical mixing option 2 (ideal)
        procedure :: reactive_mixing_iter_EI_eq_anal_ideal  !< Alternative EfI option (commented)
        procedure :: reactive_mixing_iter_EI_kin_anal_ideal_opt2       !< EI kinetics analytical option 2
        procedure :: compute_u_mix                                !< Compute tilde molalities for transport
        procedure :: compute_react_term_EE_eq_kin                   !< Explicit-Euler equilibrium-kinetics reaction iteration
        procedure :: compute_react_term_EE_kin                      !< Explicit-Euler kinetics reaction iteration
        procedure :: compute_dfk_dc_EI_ideal                          !< Compute kinetics Jacobian (EfI method)
        procedure :: compute_dfk_dc1_EI_ideal                   !< Primary kinetics Jacobian (EfI ideal)
        procedure :: Newton_EI_eq_kin_anal_ideal_opt2         !< Newton-Raphson EI equilibrium-kinetics option 2
        procedure :: Newton_EI_kin_anal_ideal_opt2            !< Newton-Raphson EI kinetics option 2
        procedure :: Newton_EI_eq_ideal                     !< General EfI kinetics Newton (commented)
    !> Update procedures for concentration and state variables
        procedure :: update_conc_aq_prim_species                    !< Update aqueous primary species concentrations
        procedure :: update_conc_prim_species                       !< Update all primary species concentrations
        procedure :: update_conc_aq_species                         !< Update all aqueous species concentrations
        procedure :: update_old_attributes               !< Update variable activity species (commented)
        procedure :: change_spec_alg_aq_chem                        !< Change speciation algorithm for aqueous chemistry
        procedure :: update_conc_nc                                 !< Update non-constant concentrations
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
        procedure :: modify_mix_ratios_reacts
    end type                                                                !< End of aqueous_chemistry_c type definition
        
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
        
        subroutine reactive_mixing_iter_EI_eq_anal_ideal(this,c1_old,c_hat,mix_ratio_r_old,mix_ratio_r_new,&
            Delta_t,theta,conc_nc,conc_comp,c1_downstream)    !< Update aqueous species concentrations (commented)
            import aqueous_chemistry_c                                   !< Import aqueous chemistry class
            implicit none                                                !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                         !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:) !> aqueous primary concentrations at previous time step
            real(kind=8), intent(in) :: c_hat(:) !> variable activity species concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r_old !> mixing ratio of reaction amount in this target
            real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of reaction amount in this target
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
            real(kind=8), intent(out) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(out) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
        end subroutine
        
        !> Component mixing iteration with ideal solution assumption
        subroutine mixing_iter_comp_ideal(this,c1_old,c_hat,mix_ratio_r_old,mix_ratio_r_new,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            !real(kind=8), intent(in) :: c2v_ig(:)                        !< Initial guess (commented - not needed for ideal)
            real(kind=8), intent(in) :: c_hat(:)                         !< Transport-advected concentrations (tilde notation)
            ! real(kind=8), intent(in) :: rk_tilde(:)                        !< Kinetic reaction rates (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r_old                    !< Mixing ratio for reactions previous time step
            real(kind=8), intent(in) :: mix_ratio_r_new                      !< Mixing ratio for reactions new time step
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
            !real(kind=8), intent(in), optional :: porosity                !< Optional porosity parameter (commented)
            !real(kind=8), intent(in), optional :: Delta_t                 !< Optional time step (commented)
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
            mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
            import aqueous_chemistry_c                                      !< Import aqueous chemistry class
            !import diag_matrix_c                                            !< Import diagonal matrix class for transport
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                            !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: c1_old(:)                          !< Primary species concentrations from previous time step
            real(kind=8), intent(in) :: c_hat(:)                          !< Transport-advected concentrations (tilde notation)
            real(kind=8), intent(in) :: mix_ratio_r                        !< Mixing ratio of kinetic reaction rate
            real(kind=8), intent(in) :: Delta_t                            !< Time step size for temporal discretization
            real(kind=8), intent(in) :: theta                              !< Time weighting factor (0=explicit, 1=implicit)
            real(kind=8), intent(inout) :: conc_nc(:)                        !< Non-constant concentrations output
            real(kind=8), intent(inout) :: conc_comp(:)                     !< Concentration components (commented)
        end subroutine
        

        subroutine reactive_mixing_iter_EI_kin_anal_ideal_opt2(this,c1_old,c_hat,&
            mix_ratio_r_old,mix_ratio_r_new,Delta_t,&
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
            real(kind=8), intent(in) :: mix_ratio_r_old !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:)
            real(kind=8), intent(inout) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            real(kind=8), intent(in), optional :: c1_downstream(:) !> primary concentrations of closest downstream target water
        end subroutine

        subroutine reactive_mixing_iter_EI_kin_anal_ideal_opt2_bis(this,c1_old,&
            c_hat,mix_ratio_r_old,mix_ratio_r_new,Delta_t,theta,conc_nc,conc_comp)
           import aqueous_chemistry_c
           import diag_matrix_c 
           !import time_discr_c
           implicit none
           class(aqueous_chemistry_c) :: this
           real(kind=8), intent(in) :: c1_old(:)
           !real(kind=8), intent(in) :: c2v_ig(:)
           real(kind=8), intent(in) :: c_hat(:)
           real(kind=8), intent(in) :: mix_ratio_r_old !> mixing ratio of kinetic reaction rate in this target
           real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction rate in this target
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

        
        subroutine compute_dfk_dc1_EI_ideal(this,c1,c2v,drk_dc,Delta_t,theta,mix_ratio_r,dfk_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c1(:)
            real(kind=8), intent(in) :: c2v(:)
            real(kind=8), intent(in) :: drk_dc(:,:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: mix_ratio_r
            real(kind=8), intent(out) :: dfk_dc1(:,:)
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
        
      
        subroutine Newton_EI_eq_kin_anal_ideal_opt2(this,u_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        
        subroutine Newton_EI_kin_anal_ideal_opt2(this,c_hat,mix_ratio_r_old,&
            mix_ratio_r_new,Delta_t,theta,conc_nc,niter,&
            CV_flag)
            import aqueous_chemistry_c
            !import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_hat(:) !> component concentrations after mixing
            real(kind=8), intent(in) :: mix_ratio_r_old !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction rate in this target
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
        
        
        subroutine compute_Re(this,cv_hat,Delta_t,lambda_r)
            import aqueous_chemistry_c
            implicit none  
            class(aqueous_chemistry_c) :: this !> object
            real(kind=8), intent(in) :: cv_hat(:) !> concentrations secondary variable activity species after mixing
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: r2nc_mix(:) !> reaction part of concentrations of secondary variable activity species after mixing
            real(kind=8), intent(in) :: lambda_r !> mixing ratio factor for reaction rates
        end subroutine

        subroutine compute_Re_kin(this,c2v_hat,Delta_t,lambda_r)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2v_hat(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio associated to "this" argument [-]
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
        
        subroutine compute_rk_new(this,rk_new)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: rk_new(:)
        end subroutine
        
        subroutine compute_Jacobian_rk_anal(this,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:,:)
        end subroutine
        
        subroutine compute_rk_Jac_rk_anal(this,rk_new,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: rk_new(:) !>  (must be allocated)
            real(kind=8), intent(out) :: drk_dc(:,:)
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
            real(kind=8), allocatable :: d_log_gamma_d_I(:) !> Derivative of log₁₀(activity coefficients) w.r.t. ionic strength: ∂(log γ)/∂I (must be already allocated) [M⁻¹]
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
    
    contains                                                                !< Beginning of implementation section
    
     
        !> Allocate memory for aqueous species concentrations
        subroutine allocate_conc_aq_species(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            allocate(this%concentrations(this%aq_phase%num_species),this%conc_old(this%aq_phase%num_species),&
                this%conc_old_old(this%aq_phase%num_species))              !< Allocate current, previous, and previous-previous concentration arrays
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
                                                                            !< Activity = concentration × activity_coefficient (a = c × γ)
        end subroutine
        
        !> Compute activities of aqueous primary species
        subroutine compute_a1(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            this%activities(this%ind_prim_species)=&
                this%concentrations(this%ind_prim_species)*&
                (10**(this%log_act_coeffs(this%ind_prim_species)))          !< Primary activities: a1 = c1 × γ1
        end subroutine
        
        !> Compute activities of secondary species
        subroutine compute_a2(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            this%activities(this%ind_sec_species)=&
                this%concentrations(this%ind_sec_species)*&
                (10**(this%log_act_coeffs(this%ind_sec_species)))           !< Secondary activities: a2 = c2 × γ2
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
        !! \note Activity coefficients must be computed beforehand (e.g., via Debye-Hückel,
        !!       Davies, or Pitzer equations). This subroutine only applies them to get activities.
        !!
        !! \warning For ideal solutions (γ = 1), activity equals concentration.
        !!
        !! \see compute_activities_diss_solids, compute_act_water, compute_ionic_strength
        !!
        subroutine compute_activities_aq(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            !> \subsection compute_dissolved_activities Compute activities of dissolved species
            !! Update activities for all dissolved species (ions and complexes) in the
            !! aqueous phase using a_i = γ_i * c_i
            call this%compute_activities_diss_solids()                     !< Apply activity coefficients to dissolved solid concentrations
            
            !> \subsection compute_water_activity Compute water activity
            !! Calculate water activity using Raoult's law: a_w ≈ 1 - 0.018 * Σm_i
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
        !!          (e.g., (Ca,Mg)CO₃), activities depend on mole fractions.
        !!
        !! \see compute_activities_aq, compute_activities_diss_solids, solid_chemistry_c::compute_activities_solids
        !!
        subroutine compute_activities(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            !> \subsection compute_aq_activities_all Compute aqueous activities
            !! Update activities for all dissolved species (ions and complexes) using
            !! the relationship a_i = γ_i * c_i
            call this%compute_activities_diss_solids()                     !< Apply activity coefficients to dissolved concentrations
            
            !> \subsection compute_water_activity_all Compute water activity
            !! Calculate water activity using Raoult's law approximation:
            !! a_w ≈ 1 - 0.018 * Σm_i for dilute solutions
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
                                                                            !< Activity calculation: a_i = c_i × γ_i
        end subroutine
        
        !> Compute water activity using Raoult's law approximation
        subroutine compute_act_water(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            if (this%aq_phase%ind_wat>0) then                              !< If water index is defined
                this%activities(this%indices_aq_species(this%aq_phase%ind_wat))=1d0-0.018*sum(this%concentrations(&
                    this%indices_aq_species(this%aq_phase%ind_diss_solids)))
                                                                            !< First-order Raoult's law: a_w = 1 - 0.018 × Σc_i
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
                !< Layout of conc_var_act_species (matches get_conc_nc / get_conc_nc_old):
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
        !! - \f$ m_i \f$ is the molality of species i [mol/kg H₂O]
        !! - \f$ z_i \f$ is the charge (valence) of species i [-]
        !! - The summation is over all aqueous species
        !!
        !! The ionic strength is a measure of the total concentration of ions in solution
        !! and is critical for calculating activity coefficients in non-ideal solutions.
        !! Higher ionic strength leads to greater deviations from ideal behavior.
        !!
        !! \note Concentrations must be in molalities (mol/kg H₂O), not molarities (mol/L).
        !!       The ionic strength is stored in this%ionic_strength.
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
            
            integer(kind=4) :: i                                            !< Loop counter for aqueous species
            
            !> \subsection init_ionic_strength Initialize ionic strength to zero
            !! Start with zero and accumulate contributions from each species.
            this%ionic_strength=0d0                                         !< Initialize ionic strength to zero [mol/kg]
            
            !> \subsection loop_species Loop through all aqueous species
            !! Calculate contribution of each species to ionic strength: c_i * z_i²
            do i=1,this%aq_phase%num_species                                !< Iterate over all aqueous species
                !> Add contribution of species i: molality × (valence)²
                !! - this%concentrations(i): molality of species i [mol/kg H₂O]
                !! - this%aq_phase%aq_species(i)%valence: charge of species i [-]
                !! - Squared valence ensures positive contribution regardless of charge sign
                this%ionic_strength=this%ionic_strength+&
                    this%concentrations(i)*this%aq_phase%aq_species(i)%valence**2
                                                                            !< Accumulate c_i * z_i² term
            end do
            
            !> \subsection scale_factor Apply the 1/2 factor
            !! The ionic strength definition includes a factor of 0.5:
            !! I = (1/2) Σ(c_i * z_i²)
            this%ionic_strength=0.5*this%ionic_strength                     !< Multiply by 1/2 to get final ionic strength [mol/kg]
        end subroutine
        
        !> \brief Computes logarithmic equilibrium constants for all equilibrium reactions
        !!
        !! \details This subroutine updates the log₁₀(K) values for all equilibrium reactions
        !! in the chemical system based on the current temperature. The equilibrium constant
        !! K is temperature-dependent according to the van't Hoff equation:
        !! \f[
        !! \frac{d \ln K}{dT} = \frac{\Delta H°}{RT^2}
        !! \f]
        !! where:
        !! - K is the equilibrium constant [-]
        !! - T is absolute temperature [K]
        !! - ΔH° is the standard enthalpy of reaction [J/mol]
        !! - R is the universal gas constant [J/(mol·K)]
        !!
        !! **Algorithm**:
        !! 1. Update logK for all equilibrium reactions in the chemical system
        !! 2. If solid chemistry is present, update logK for reactions in reactive zone
        !!
        !! \note The temperature is taken from this%temp which must be set beforehand.
        !!       The computed logK values are stored within each eq_reaction object.
        !!
        !! \warning This subroutine modifies the equilibrium constants of the chemical system.
        !!          It must be called whenever temperature changes before performing speciation.
        !!
        !! \see eq_reaction_c::compute_logK_dep_T, set_temp
        !!
        subroutine compute_log_K_aq_chem(this)
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< [in,out] Aqueous chemistry object
            
            integer(kind=4) :: i                                            !< Loop counter for equilibrium reactions
            
            !> \subsection update_all_reactions Update logK for all equilibrium reactions
            !! Iterate through all equilibrium reactions in the chemical system and
            !! recompute their equilibrium constants at the current temperature.
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_eq_reactions !< Loop over all equilibrium reactions
                !> Compute temperature-dependent log₁₀(K) for reaction i
                !! Uses thermodynamic data (ΔH°, ΔS°) stored in each reaction object
                !! to calculate K(T) from the van't Hoff equation
                call this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(i)%compute_logK_dep_T(this%temp)
                                                                            !< Update logK for reaction i at current temperature
            end do
            
            !> \subsection update_reactive_zone Update logK for reactive zone reactions (if applicable)
            !! If solid chemistry is associated (minerals present), also update the
            !! equilibrium constants for reactions specific to the reactive zone.
            !! This handles heterogeneous reactions (mineral dissolution/precipitation).
            if (associated(this%solid_chemistry)) then                     !< Check if solid chemistry exists
                !> Loop through reactive zone equilibrium reactions
                !! The reactive zone may have a subset of the total chemical system reactions
                !! that are relevant for the current spatial location or time step
                do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !< Reactive zone eq reactions
                    !> Get reaction index and update its logK
                    !! ind_eq_reacts(i) maps from reactive zone index to chemical system index
                    call this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(&
                        this%solid_chemistry%reactive_zone%ind_eq_reacts(i))%compute_logK_dep_T(this%temp)
                                                                            !< Update logK for reactive zone reaction i
                end do
            end if
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
            
            !> \subsection validate_reactive_zone Validate solid chemistry initialization
            !! Check that the solid chemistry object has been properly initialized with
            !! a reactive zone before creating the association. The reactive zone contains
            !! information about which minerals are present and their reaction kinetics.
            if (associated(solid_chemistry%reactive_zone)) then            !< Check if reactive zone is associated
                !> Create pointer association to solid chemistry
                !! This allows the aqueous chemistry object to access mineral data,
                !! equilibrium reactions involving minerals, and kinetic rate laws
                this%solid_chemistry=>solid_chemistry                      !< Point to the solid chemistry object
            else
                !> Terminate if reactive zone is not associated
                !! This is a critical error because mineral reactions cannot be computed
                !! without reactive zone information
                error stop "Solid chemistry object is not associated to a reactive zone"
                                                                            !< Fatal error: invalid solid chemistry
            end if
        end subroutine

        !> \brief Associates a solid chemistry object from the previous time step
        !!
        !! \details Creates a pointer association to a solid chemistry object
        !! representing the state at the previous time step. If not provided,
        !! defaults to the current solid chemistry.
        !!
        !! @param solid_chemistry_old optional solid chemistry object from previous time step;
        !!        defaults to current solid chemistry
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
            end if
        end subroutine
        
        !> \brief Associates a gas chemistry object to this aqueous chemistry
        !!
        !! \details This subroutine creates a pointer association between the aqueous
        !! chemistry object and a gas chemistry object (containing gases and their
        !! reactive zones). This coupling is necessary for:
        !! - Gas dissolution/exsolution reactions (e.g., CO₂, O₂, CH₄)
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
        !! \details Allocates Rk, rk_old, rk_new, rk_mean, Rk_est, Rk_accum,
        !! re_mean, and Re arrays for aqueous, solid, and gas chemistry objects.
        !! All arrays are initialised to zero.
        !!
        subroutine allocate_reaction_rates(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (.not. associated(this%solid_chemistry)) then
                error stop "solid_chemistry pointer not associated in allocate_reaction_rates"
            end if
            if (.not. associated(this%solid_chemistry%reactive_zone)) then
                error stop "reactive_zone pointer not associated in allocate_reaction_rates"
            end if
            if (associated(this%solid_chemistry%reactive_zone%chem_syst)) then
                associate(n_kin => this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                                   this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts)
                allocate(this%Rk(n_kin))
                allocate(this%rk_old(n_kin))
                allocate(this%rk_old_old(n_kin))
                allocate(this%rk_old_old_old(n_kin))
                allocate(this%rk_mean(n_kin))
                allocate(this%Rk_est(n_kin))
                allocate(this%rk_new(n_kin))
                allocate(this%Rk_accum(n_kin))
                end associate
                this%Rk=0d0 !> by default
                this%Rk_est=0d0 !> by default
                this%rk_old=0d0 !> by default
                this%rk_new=0d0 !> by default
                this%rk_old_old=0d0 !> by default
                this%rk_old_old_old=0d0 !> by default
                this%rk_mean=0d0 !> by default
                this%Rk_accum=0d0 !> by default
            else
                error stop "Chemical system not associated to reactive zone"
            end if
            allocate(this%re_mean(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%Re(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            this%re_mean=0d0 !> by default
            this%Re=0d0 !> by default
            if (associated(this%solid_chemistry)) then
                if (.not. allocated(this%solid_chemistry%re_mean)) then
                    allocate(this%solid_chemistry%re_mean(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    this%solid_chemistry%re_mean=0d0 !> by default
                    this%solid_chemistry%Re=0d0 !> by default
                end if
                if (.not. allocated(this%solid_chemistry%Rk)) then
                    associate(n_min_kin => this%solid_chemistry%mineral_zone%num_minerals_kin)
                    allocate(this%solid_chemistry%Rk(n_min_kin))
                    allocate(this%solid_chemistry%rk_old(n_min_kin))
                    allocate(this%solid_chemistry%rk_old_old(n_min_kin))
                    allocate(this%solid_chemistry%rk_old_old_old(n_min_kin))
                    allocate(this%solid_chemistry%rk_mean(n_min_kin))
                    allocate(this%solid_chemistry%rk_new(n_min_kin))
                    allocate(this%solid_chemistry%Rk_est(n_min_kin))
                    allocate(this%solid_chemistry%Rk_accum(n_min_kin))
                    end associate
                    this%solid_chemistry%rk_new=0d0 !> by default
                    this%solid_chemistry%Rk_est=0d0 !> by default
                    this%solid_chemistry%rk_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old_old=0d0 !> by default
                    this%solid_chemistry%rk_mean=0d0 !> by default
                    this%solid_chemistry%Rk=0d0 !> by default
                    this%solid_chemistry%Rk_accum=0d0 !> by default
                end if
            else
                error stop "Solid chemistry not associated to aqueous chemistry"
            end if
            if (associated(this%gas_chemistry)) then
                if (.not. allocated(this%gas_chemistry%re_mean)) then
                    allocate(this%gas_chemistry%re_mean(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%Re(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    this%gas_chemistry%re_mean=0d0 !> by default
                    this%gas_chemistry%Re=0d0 !> by default
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
            log_a=(log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%eq_cst)-dot_prod)/this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(ind_eq_rct)%stoichiometry(ind_act_react)
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
        !! Mathematical Context: Salinity = TDS/(1+TDS) where TDS = Σ(ci × MWi) [kg_solute/kg_water]
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
            !! mass contributions: mass_i = molality_i × molecular_weight_i
            do i=1,this%aq_phase%num_species-this%aq_phase%wat_flag         !< Loop over dissolved solids (water excluded)
                !> Add mass contribution for species i
                !! Extract concentration from full array using index mapping,
                !! then multiply by molecular weight to convert mol/kg → kg/kg
                TDS=TDS+this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids(i))) * &
                    this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i))%molecular_weight
                                                                            !< TDS += m_i × MW_i [kg_solute/kg_water]
            end do
            
            !> \subsection compute_total_mass Calculate total fluid mass
            !! Total mass = mass of water (1 kg basis) + mass of all dissolved solutes
            fluid_mass=1d0+TDS                                             !< Total fluid mass = 1 kg water + TDS [kg]
            
            !> \subsection compute_salinity_value Compute salinity as mass fraction
            !! Salinity is the dimensionless ratio of solute mass to total fluid mass
            this%salinity=TDS/fluid_mass                                   !< S = TDS/(1+TDS) [-]
        end subroutine
        
        !> Compute alkalinity as charge balance contribution from all aqueous species
        !! Mathematical Context: Alkalinity = Σ(ci × alk_contrib_i) where alk_contrib is species charge contribution
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
                                                                            !< Add contribution: concentration × alkalinity_factor
            end do
            this%alkalinity = alkalinity                                   !< Store computed alkalinity in object
        end subroutine
        
        !> Convert molarities to molalities using density and salinity corrections
        !! Mathematical Context: molality = molarity / (density × (1-salinity))
        !! This accounts for volume changes due to dissolved solids and temperature effects
        subroutine compute_molalities(this)                                 !< Convert from molarities to molalities
            implicit none                                                   !< No implicit variable declarations
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            
            ! real(kind=8) :: log_a,dot_prod                                !< Commented legacy variables
            ! integer(kind=4) :: ind_act_react,i                           !< Commented legacy variables
            
            !> Convert dissolved solids concentrations from molarity to molality
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))/(this%density*(1d0-this%salinity))
                                                                            !< molality = molarity / (density × (1-salinity))
            !> Convert water concentration using density correction only
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))/this%density
                                                                            !< water_molality = water_molarity / density
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
        function get_c1_old_old(this) result(c1_old_old)                                   !< Function returns primary concentrations array
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
        function get_conc_nc(this) result(conc_nc) !> gets conc var act species per unit volume of water
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

        !> \brief Gets concentrations of variable activity species from previous time step
        !!
        !! \details Assembles and returns variable activity species concentrations
        !! from the previous time step. Same ordering as get_conc_nc.
        !!
        !! @return conc_nc variable activity species concentrations from previous time step
        !!
        function get_conc_nc_old(this) result(conc_nc) !> gets conc var act species per unit volume of water from previous time step
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
                conc_nc(i)=this%conc_old(this%ind_var_act_species(i))
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
                !! this%conc_old(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=this%conc_old(&
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                        this%solid_chemistry%conc_old(this%solid_chemistry%mineral_zone%num_minerals+1)
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=this%solid_chemistry%conc_old(&
                        this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%num_solids)
                end if
            if (associated(this%gas_chemistry)) then
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act)=this%gas_chemistry%conc_old(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
            end if
        end function
        
        !> \brief Get logarithmic activity coefficients of secondary variable activity species.
        !>
        !> This function retrieves the logarithmic activity coefficients (log₁₀ γ) for all secondary
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
        !> \return log_gamma2v Vector of log₁₀ γ₂ₙ꜀ values [-] (size = num_eq_reactions)
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
        !> This function retrieves the logarithmic activity coefficients (log₁₀ γ) for all species
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
        !> 9. Equilibrium constant activity minerals (log γ = 0, pure solid phases)
        !> 10. Surface complexes (adsorbed species on mineral surfaces)
        !> 11. Equilibrium variable activity gases (ideal gas behavior)
        !> 12. Equilibrium constant activity gases (ideal gas behavior)
        !> This comprehensive ordering allows reconstruction of the complete thermodynamic state
        !> for reactive transport calculations.
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state and activity coefficients
        !> \return log_gamma_vec Vector of log₁₀ γ values [-] for all species (size = num_species)
        !>
        !> \note For constant activity species (pure phases, ideal gases), log γ = 0 by definition.
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
            
            log_res=log10(this%solid_chemistry%reactive_zone%get_eq_csts_react_zone())-&
            matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,1:this%aq_phase%num_species),log10(this%activities))
            if (inf_norm_vec_real(log_res)>this%solid_chemistry%reactive_zone%CV_params%log_rel_tol) then
                print *, inf_norm_vec_real(log_res)
                error stop "Error in aqueous activities"
            end if
        end subroutine
        
        !> \brief Gets aqueous variable activity species concentrations after applying a permutation
        !!
        !! @param perm permutation array
        !! @return c_nc_aq permuted aqueous concentrations
        !!
        function get_conc_nc_aq_perm(this,perm) result(c_nc_aq) !> gets concentrations of aqueous species after a permutation
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
        
        function compute_saturation_kin_min(this,ind_min) result(saturation)
        !> This function is meant to be used for minerals NOT in equilibrium
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: ind_min !> index of kinetic mineral in chemical system
            real(kind=8) :: saturation !> saturation index
            
            real(kind=8) :: IAP
            integer(kind=4) :: i,ind_kin_react

            IAP=1d0
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(ind_min)%num_aq_rk
                IAP=IAP*this%activities(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                    ind_min)%indices_aq_phase(i)))**(-&
                    this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(ind_min)%stoichiometry(&
                    this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(ind_min)%indices_react_species(i)))
            end do
            saturation=IAP*this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(ind_min)%eq_cst
        end function
        
     
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
        !! Mathematical Context: I = 1/2 * Σ(ci * zi²) where I is ionic strength, ci is concentration, zi is charge
        !! The derivative dI/dci = 1/2 * zi² + Σ(dcsec/dci * zsec²) for secondary species dependencies
        subroutine compute_dI_dc1(this,dc2aq_dc1,dI_dc1)
            class(aqueous_chemistry_c) :: this                             !< Aqueous chemistry object (polymorphic)
            real(kind=8), intent(in) :: dc2aq_dc1(:,:)                     !< Jacobian of secondary aqueous concentrations w.r.t primary
            real(kind=8), intent(out) :: dI_dc1(:)                         !< Jacobian of ionic strength w.r.t primary concentrations
                                                                            !! (must be already allocated)
            
            integer(kind=4) :: i                                           !< Loop index for primary species
            
            !> For each primary species, compute dI/dci = 1/2*zi² + Σ(dcsec/dci * zsec²)
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                dI_dc1(i)=5d-1*this%aq_phase%aq_species(i)%valence**2 + &  !< Direct contribution: 1/2 * zi²
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
                !> Compute Newton residual Jacobian with mixing ratio: df/dc = I - θ*Δt*mix_ratio*U*S*dR/dc
                dfk_dc=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)-theta*Delta_t*mix_ratio_Rk*&
                    matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc)  !< Apply stoichiometric transformation with mixing
            end subroutine

        !> Compute aqueous component concentrations after mixing (transport step)
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
             print *, "Warning: concentration of aqueous primary species is zero"
            print *, "Aqueous primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_prim_species(i)))%name
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
            print *, "Warning: concentration of primary species is zero"
            print *, "Primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(i)))%name
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
    c_mix = this%get_conc_nc_old()
    num_perm = 0
    do i = 1, n_aq_vas
        idx_aq = ind_aq(ind_diss(i))
        !> First water contribution from "this"
        c_mix_aq(i) = mixing_ratios_conc(1) * this%conc_old(idx_aq)
        !> Build permutation table where ind_var_act_species differs from indices_aq_species
        if (ind_vas(i) /= idx_aq) then
            !> Check if i is already recorded in perm(:,2) — inline search
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

!> Compute aqueous mixed concentrations using pre-cached global contiguous arrays.
!! all_conc_old(:, w) stores conc_old of water w in that water's own species order.
!! all_ind_aq_sp(k, w) stores indices_aq_species(k) of water w,
!! i.e. the position in water w's conc_old of aq_phase species k.
!! all_conc_old and all_ind_aq_sp do NOT contain data for "this" (only mixing waters).
!! Output c_mix contains all aqueous species after mixing, ordered as:
!!   1. Primary species (1..n_aq_prim)
!!   2. Secondary variable activity species
!!   3. Secondary constant activity species
!! following the reactive zone associated to "this".
subroutine compute_c_mix_global(this, all_conc_old, all_ind_aq_sp, mix_col, &
    num_mix_wat, mixing_ratios_conc, c_mix)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: this !> target water (structural reference)
    real(kind=8), intent(in) :: all_conc_old(:,:) !> (max_n_conc, num_waters) conc_old cache for mixing waters (excludes "this")
    integer(kind=4), intent(in) :: all_ind_aq_sp(:,:) !> (num_aq_species, num_waters) indices_aq_species for mixing waters (excludes "this")
    integer(kind=4), intent(in) :: mix_col(:) !> mixing water indices (1..num_mix_wat), excludes "this"
    integer(kind=4), intent(in) :: num_mix_wat !> number of mixing waters (upstream + downstream, excludes "this")
    real(kind=8), intent(in) :: mixing_ratios_conc(:) !> mixing ratios (1=self, 2..num_mix+1=mixing waters)
    real(kind=8), intent(out), allocatable :: c_mix(:) !> mixed aqueous species concentrations (primary, sec var act, sec cst act)
!> Variables
    integer(kind=4) :: i, j, w, aq_idx
    integer(kind=4) :: n_aq, n_aq_prim, n_sec, n_vas, n_prim, n_aq_sec
    
    n_aq = this%aq_phase%num_species
    n_aq_prim = size(this%ind_prim_species)
    n_sec = size(this%ind_sec_species)
    !> Full variable-activity vector layout (matches get_conc_nc_old):
    !>   1..n_aq_prim                   : aqueous primary
    !>   n_aq_prim+1..n_prim            : non-aqueous primary (e.g. free surface site)
    !>   n_prim+1..n_prim+n_aq_sec      : aqueous secondary variable activity
    !>   n_prim+n_aq_sec+1..n_vas       : non-aqueous secondary var act (e.g. surface complexes)
    n_vas    = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    n_prim   = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_aq_sec = this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
    
    block
    real(kind=8) :: c_tmp(n_aq)
    !> Step 1: Mix in the order of "this" concentrations attribute
    !> Self contribution comes from this%conc_old directly
    !> indices_aq_phase(i) gives the aq_phase index of species at position i in concentrations
    !> all_ind_aq_sp(aq_idx, w) = position in water w's conc_old of aq_phase species aq_idx
    do i = 1, n_aq
        aq_idx = this%indices_aq_phase(i)
        c_tmp(i) = mixing_ratios_conc(1) * this%conc_old(i)
        do j = 1, num_mix_wat
            w = mix_col(j)
            c_tmp(i) = c_tmp(i) + mixing_ratios_conc(j+1) * all_conc_old(all_ind_aq_sp(aq_idx, w), w)
        end do
    end do
    !> Step 2: Allocate c_mix over the full variable-activity vector and
    !>         prefill with the previous step's concentrations so that the
    !>         non-aqueous slots (surface complexes / minerals / gases) carry
    !>         their old values (mixing only acts on aqueous species).
    allocate(c_mix(n_vas))
    c_mix = this%get_conc_nc_old()
    !> Step 3: Overwrite aqueous slots from the mixed values, reordered from
    !>         "concentrations" order to (aq primary, aq secondary) order.
    do i = 1, n_aq_prim
        c_mix(i) = c_tmp(this%ind_prim_species(i))
    end do
    do i = 1, n_aq_sec
        c_mix(n_prim + i) = c_tmp(this%ind_sec_species(i))
    end do
    end block
end subroutine compute_c_mix_global

!> This subroutine computes the concentration of aqueous variable activity species & aqueous kinetic reaction rates after water mixing
!! We assume all kinetic reactions are aqueous
!! we assume all mixing waters have the same kinetic reactions
subroutine compute_c_and_rk_tilde(mix_waters,mixing_ratios_conc,mixing_ratios_R,c_mix,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
    real(kind=8), intent(in) :: mixing_ratios_R(:)  !> same order as "mix_waters"
    real(kind=8), intent(out) :: c_mix(:) !> concentration of aqueous variable activity species after mixing
    real(kind=8), intent(out) :: rk_tilde(:) !> aqueous kinetic reaction rates after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: c_mix_aux(:) !> auxiliary variable to store the concentration of aqueous variable activity species after mixing
    logical :: flag !> flag to check if an index is in "perm"
!> Pre-process
    allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
!> Process
    c_mix=0d0 !> initialise concentration of aqueous variable activity species after mixing
    rk_tilde=0d0 !> initialise aqueous kinetic reaction rates after mixing
    num_perm=0 !< counter in "perm"
    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
        !> Autentica chapuza
        if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
            if (num_perm>0) then
                call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
                if (flag .eqv. .true.) then
                    continue
                else
                    num_perm=num_perm+1
                    perm(num_perm,1)=i
                    perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                        mix_waters(1)%ind_var_act_species)
                end if
            else
                num_perm=num_perm+1
                perm(num_perm,1)=i
                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                    mix_waters(1)%ind_var_act_species)
            end if
        end if
        do j=1,num_mix_wat
            c_mix(i)=c_mix(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
                mix_waters(j)%aq_phase%ind_diss_solids(i)))
            
        end do
    end do
    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        do j=1,num_mix_wat
            rk_tilde(i)=rk_tilde(i)+mixing_ratios_R(j)*mix_waters(j)%rk(i)
        end do
    end do
    if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
            j=1
            k=1
            i=1
            do
                if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
                    mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
                    rk_tilde(i)=rk_tilde(i)+mixing_ratios_R(j)*mix_waters(j)%solid_chemistry%rk(i)
                    if (j<num_mix_wat) then
                        j=j+1
                        i=1
                        k=1
                    else
                        exit
                    end if
                else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
                    k=k+1
                else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
                    i=i+1
                    k=1
                else if (j<num_mix_wat) then
                    j=j+1
                    i=1
                    k=1
                else
                    exit
                end if
            end do
    end if
    c_mix_aux=c_mix !> we store the concentration of aqueous variable activity species after mixing (chapuza)
    do i=1,num_perm
        !> we permute the indices in "perm"
        c_mix(perm(i,1))=c_mix_aux(perm(i,2))
        c_mix(perm(i,2))=c_mix_aux(perm(i,1))
    end do
end subroutine 

!> This subroutine computes the reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters, for the first time step
!> The estimation of the downstream waters reaction rates is done by taking the reaction rate of the closest upstream water
subroutine compute_r_tilde_impl_opt1(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
                                     mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> target waters that mix with each other, ordered in upstream & downstream waters
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> same order as "mix_waters"
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: Sk_nc_up(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_up_old(:,:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_old(:,:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: Se_nc(:,:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: U_SeT(:,:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_new(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_up_old(:) !> auxiliary variable to store the reaction rate contributions after mixing
    real(kind=8) :: tol !> tolerance
    
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt1"
    
!> Process
    tol=mix_waters(1)%solid_chemistry%reactive_zone%CV_params%zero
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !    !> Autentica chapuza
    !        
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                Sk_nc_old=mix_waters(j)%get_Sk_nc_old() !> we get the kinetic stoichiometric matrix in (j-1)-th upstream water (chapuza)
                rk_new=mix_waters(j)%get_rk_new()
                rk_old=mix_waters(j)%get_rk_old()
                r_tilde=r_tilde+mixing_ratios_reacts(j)*(theta*matmul(&
                    transpose(Sk_nc),rk_new)+(1d0-theta)*matmul(transpose(&
                    Sk_nc_old),rk_old))
            end if
            if (mix_waters(j)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
                Se_nc=mix_waters(j)%solid_chemistry%reactive_zone%get_Se_nc_react_zone(ind_var_act_species_rz) !> we get the equilibrium stoichiometric matrix in reactive zone
                allocate(U_SeT(size(comp_mat_rz,1),size(Se_nc,1)))
                U_SeT=matmul(comp_mat_rz,transpose(Se_nc)) !> we compute U*Se_nc^T
                if (norm_mat_inf(U_SeT)>=tol) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Se_nc),&
                        mix_waters(j)%get_re_mean())
                end if
                deallocate(U_SeT)
            end if
        end do
    !> Downstream waters
        if (num_up>0) then
            Sk_nc_up=Sk_nc !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
            Sk_nc_up_old=Sk_nc_old !> we get the kinetic stoichiometric matrix in (j-1)-th upstream water (chapuza)
            rk_up_old=rk_old !> we get the kinetic stoichiometric matrix in (j-1)-th upstream water (chapuza)
            do j=1,num_down
                if (mix_waters(num_up+j)%indices_rk%num_cols>0) then
                    Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    Sk_nc_old=mix_waters(num_up+j)%get_Sk_nc_old() !> we get the kinetic stoichiometric matrix in (j-1)-th upstream water (chapuza)
                    if (norm_mat_inf(Sk_nc-Sk_nc_up)<tol) then
                        r_tilde=r_tilde+mixing_ratios_reacts(num_up+j-1)*(theta*matmul(&
                            transpose(Sk_nc_up),rk_new)+(1d0-theta)*matmul(transpose(&
                            Sk_nc_up_old),rk_up_old))
                    else
                        rk_old=mix_waters(num_up+j)%get_rk_old()
                        r_tilde=r_tilde+mixing_ratios_reacts(num_up+j)*matmul(transpose(Sk_nc_old),rk_old)
                    end if
                end if
            end do
        else
            do j=1,num_down
                if (mix_waters(j)%indices_rk%num_cols>0) then
                    Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    Sk_nc_old=mix_waters(j)%get_Sk_nc_old() !> we get the kinetic stoichiometric matrix in (j-1)-th upstream water (chapuza)
                        rk_old=mix_waters(j)%get_rk_old()
                        r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Sk_nc_old),rk_old)
                end if
            end do
        end if
    ! !> Downstream waters
        ! !> Current water
end subroutine

!> This subroutine computes the kinetic reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by extrapolation with the two previous reaction rates
subroutine compute_r_tilde_impl_opt2(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> indices of variable activity species in the current reactive zone
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
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: rk_est(:) !> kinetic stoichiometric matrix associated to the previous target
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    !real(kind=8) :: mu=1d0 !> parameter for extrapolation of downstream waters rk
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down+1) error stop "Dimension error in compute_r_tilde_impl_opt2"
    
!> Process
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    !rk_tilde_down=0d0 !> initialise kinetic reaction contributions after mixing
    r_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j+1)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_R(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Sk_nc),theta*mix_waters(j+1)%rk_new+&
                        (1d0-theta)*mix_waters(j+1)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j+1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Sk_nc),theta*mix_waters(j+1)%solid_chemistry%rk_new+&
                        (1d0-theta)*mix_waters(j+1)%solid_chemistry%rk_old)
                end if
                deallocate(Sk_nc)
            end if
        end do
    !> Downstream waters
        do j=1,num_down
            if (mix_waters(num_up+j+1)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                if (mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    allocate(rk_est(&
                        mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts))
                    rk_est=(1d0+&
                        mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%CV_params%est_prm*theta)*&
                        mix_waters(num_up+j+1)%rk_old-&
                        mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%CV_params%est_prm*&
                        theta*mix_waters(num_up+j+1)%rk_old_old
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j+1)*&
                        matmul(transpose(Sk_nc),rk_est)
                    call mix_waters(num_up+j+1)%set_Rk_est(rk_est*Delta_t)
                    !end if
                    deallocate(rk_est)
                end if
                if (allocated(mix_waters(num_up+j+1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    allocate(rk_est(&
                        mix_waters(num_up+j+1)%solid_chemistry%mineral_zone%num_minerals_kin))
                    rk_est=(1d0+&
                        mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%CV_params%est_prm*theta)*&
                        mix_waters(num_up+j+1)%solid_chemistry%rk_old-&
                        mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%CV_params%est_prm*&
                        theta*mix_waters(num_up+j+1)%solid_chemistry%rk_old_old
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j+1)*&
                        matmul(transpose(Sk_nc),rk_est)
                    call mix_waters(num_up+j+1)%solid_chemistry%set_Rk_est(rk_est*Delta_t)
                    !end if
                    deallocate(rk_est)
                end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
        !> Current water
        if (mix_waters(1)%indices_rk%num_cols>0) then
            Sk_nc=mix_waters(1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in 1-th upstream water (chapuza)
            !rk_tilde=rk_tilde+mixing_ratios_R(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
            !deallocate(Sk_nc)
            if (mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc),&
                    mix_waters(1)%rk_old)
                !deallocate(Sk_nc)
            end if
            if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc),&
                    mix_waters(1)%solid_chemistry%rk_old)
            end if
            !Sk_nc_old=Sk_nc !> chapuza
            deallocate(Sk_nc)
        end if
end subroutine 

!> This subroutine computes the reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by assuming there is a fraction of the reaction that occurs at the original target
subroutine compute_r_tilde_impl_opt4(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> target waters that mix with each other (must be ordered in upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> component matrix in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    !integer(kind=4), intent(in) :: Lagr_flag !> 0: Eulerian, 1: Lagrangian
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_v !> number of variable activity species in downstream water
    integer(kind=4) :: n_aq_kin !> number of aqueous kinetic reactions in downstream water
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
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_r_tilde_impl_opt4"
    
!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            S_nc=mix_waters(j)%get_S_nc(ind_var_act_species_rz) !> we get the new stoichiometric matrix in j-th upstream water (chapuza)
            S_nc_old=mix_waters(j)%get_S_nc_old(ind_var_act_species_rz) !> we get the old stoichiometric matrix in (j)-th upstream water (chapuza)
            r_new=mix_waters(j)%get_r_new() !> we get the new reaction rates
            r_old=mix_waters(j)%get_r_old() !> we get the old reaction rates
            r_tilde=r_tilde+mixing_ratios_reacts(j)*(theta*matmul(&
                transpose(S_nc),r_new)+(1d0-theta)*matmul(transpose(S_nc_old),&
                r_old))
        end do
        !> Downstream waters (we estimate the reaction rates by assuming a fraction of the reaction occurs at the original target)
        do j=1,num_down
            S_nc_old=mix_waters(j+num_up)%get_S_nc_old(ind_var_act_species_rz) !> we get the old stoichiometric matrix in (j)-th upstream water (chapuza)
            r_old=mix_waters(j+num_up)%get_r_old() !> we get the old reaction rates
            r_tilde=r_tilde+mixing_ratios_reacts(j+num_up)*(1d0-theta)*matmul(transpose(S_nc_old),r_old)
                !     !> Esto de abajo esta mal (todavia no has calculado mix_waters(1)%rk_new)
                !     !> Esto de abajo esta mal (todavia no has calculado mix_waters(1)%solid_chemistry%rk_new)
        end do
        ! !> Current water
end subroutine

!> \brief Computes kinetic reaction rate contributions after water mixing implicitly (option 4 variant)
!!
!! \details Assembles the reaction rate tilde vector for mixed waters using
!! implicit time weighting with stoichiometric matrices from current and
!! upstream/downstream waters.
!!
!! @param mix_waters array of waters that mix (ordered: current, upstream, downstream)
!! @param ind_var_act_species_rz indices of variable activity species in the reactive zone
!! @param comp_mat_rz component matrix in the current reactive zone
!! @param mixing_ratios_reacts mixing ratios (same order as mix_waters)
!! @param num_up number of upstream waters
!! @param num_down number of downstream waters
!! @param theta time weighting factor
!! @param Delta_t time step
!! @param r_tilde reaction rate contributions after mixing (output)
!!
subroutine compute_r_tilde_impl_opt4_bis(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other (must be ordered in current, upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> component matrix in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> ordered in current, upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    !integer(kind=4), intent(in) :: Lagr_flag !> 0: Eulerian, 1: Lagrangian
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: ind_up !> index of upstream water
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_v !> number of variable activity species in downstream water
    integer(kind=4) :: n_aq_kin !> number of aqueous kinetic reactions in downstream water
    integer(kind=4) :: n_aq_kin_up !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin_up !> number of aqueous kinetic reactions in current water
    real(kind=8), allocatable :: rk_est(:) !> mean kinetic reaction amounts
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_up(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: A(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: b(:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: conc_nc(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: conc_nc_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: conc_nc_est(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down+1) error stop "Dimension error in compute_r_tilde_impl_opt4"
    
!> Process
    r_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up-1
            if (mix_waters(j+2)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j+2)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                Sk_nc_up=mix_waters(j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                n_aq_kin=mix_waters(j+2)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
                n_min_kin=mix_waters(j+2)%solid_chemistry%mineral_zone%num_minerals_kin
                n_min_kin_up=mix_waters(j+1)%solid_chemistry%mineral_zone%num_minerals_kin
                if (n_aq_kin>0) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Sk_nc(1:n_aq_kin,:)),&
                        theta*mix_waters(j+2)%rk_new+&
                        (1d0-theta)*mix_waters(j+2)%rk_old)
                end if
                if (allocated(mix_waters(j+2)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*(matmul(transpose(Sk_nc(&
                        n_aq_kin+1:n_aq_kin+n_min_kin,:)),&
                        theta*mix_waters(j+2)%solid_chemistry%rk_new)+matmul(&
                        transpose(Sk_nc_up(n_aq_kin+1:n_aq_kin+n_min_kin_up,:)),&
                        (1d0-theta)*mix_waters(j+1)%solid_chemistry%rk_old))
                end if
                deallocate(Sk_nc)
            end if
        end do
        !> Downstream waters (we estimate the reaction rates by assuming a fraction of the reaction occurs at the original target)
        do j=1,num_down-1
            if (j==1) then
                ind_up=1
            else
                ind_up=num_up+j
            end if
            if (mix_waters(ind_up)%indices_rk%num_cols>0) then
                Sk_nc_up=mix_waters(ind_up)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                Sk_nc=mix_waters(1+j+num_up)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in current water
                n_v=mix_waters(ind_up)%solid_chemistry%reactive_zone%chem_syst%num_var_act_species
                n_aq_kin=mix_waters(1+j+num_up)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
                n_min_kin=mix_waters(ind_up)%solid_chemistry%mineral_zone%num_minerals_kin
                if (n_aq_kin>0) then
                    r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(j+num_up)*matmul(transpose(&
                        Sk_nc(1:n_aq_kin,:)),&
                        mix_waters(1+j+num_up)%rk_old)
                    !> Esto de abajo esta mal (todavia no has calculado mix_waters(1)%rk_new)
                end if
                if (allocated(mix_waters(ind_up)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(j+num_up)*matmul(transpose(&
                        Sk_nc_up(n_aq_kin+1:n_aq_kin+n_min_kin,:)),&
                        mix_waters(ind_up)%solid_chemistry%rk_old)
                    !> Esto de abajo esta mal (todavia no has calculado mix_waters(1)%solid_chemistry%rk_new)
                end if
                deallocate(Sk_nc)
                deallocate(Sk_nc_up)
            end if
        end do
        !> Current water
        if (mix_waters(1)%indices_rk%num_cols>0) then
            Sk_nc=mix_waters(1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
            Sk_nc_up=mix_waters(1+num_up)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
            if (mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc),&
                    mix_waters(1)%rk_old)
            end if
            if (allocated(mix_waters(1+num_up)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc_up),&
                    mix_waters(1+num_up)%solid_chemistry%rk_old)
            end if
            deallocate(Sk_nc)
            deallocate(Sk_nc_up)
        end if
end subroutine

!> \brief Computes kinetic reaction rate contributions after water mixing implicitly (option 1 variant)
!!
!! \details Assembles the reaction rate tilde vector for mixed waters using
!! implicit time weighting with stoichiometric matrices. Uses a different
!! assembly strategy from option 4.
!!
!! @param mix_waters array of target waters that mix (ordered: current, upstream, downstream)
!! @param ind_var_act_species_rz indices of variable activity species in the reactive zone
!! @param comp_mat_rz component matrix in the current reactive zone
!! @param mixing_ratios_reacts mixing ratios (same order as mix_waters)
!! @param num_up number of upstream target waters
!! @param num_down number of downstream target waters
!! @param theta time weighting factor
!! @param Delta_t time step
!! @param r_tilde reaction rate contributions after mixing (output)
!!
subroutine compute_r_tilde_impl_opt1_bis(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Arguments
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> target waters that mix with each other (must be ordered in current, upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> component matrix in the current reactive zone
    real(kind=8), intent(in) :: mixing_ratios_reacts(:)  !> same order as "mix_waters"
    integer(kind=4), intent(in) :: num_up  !> number of upstream target waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream target waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    !integer(kind=4), intent(in) :: Lagr_flag !> 0: Eulerian, 1: Lagrangian
    real(kind=8), intent(out) :: r_tilde(:) !> reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: ind_up !> index of upstream water
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4) :: n_v !> number of variable activity species in downstream water
    integer(kind=4) :: n_aq_kin !> number of aqueous kinetic reactions in downstream water
    integer(kind=4) :: n_aq_kin_up !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin !> number of aqueous kinetic reactions in current water
    integer(kind=4) :: n_min_kin_up !> number of aqueous kinetic reactions in current water
    real(kind=8), allocatable :: rk_est(:) !> mean kinetic reaction amounts
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_up(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_old(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_up_old(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: rk_new(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: rk_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_up_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: U_SeT(:,:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), allocatable :: Se_nc(:,:) !> equilibrium stoichiometric matrix associated to each target
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    !real(kind=8) :: mu=0d0 !> parameter for extrapolation of downstream waters rk
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down+1) error stop "Dimension error in compute_r_tilde_impl_opt1_bis"
    
!> Process
    r_tilde=0d0 !> initialise reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j+1)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                Sk_nc_old=mix_waters(j+1)%get_Sk_nc_old() !> chapuza
                rk_new=mix_waters(j+1)%get_rk_new() !> chapuza
                rk_old=mix_waters(j+1)%get_rk_old() !> chapuza
                r_tilde=r_tilde+mixing_ratios_reacts(j+1)*(&
                    theta*(matmul(transpose(Sk_nc),theta*rk_new))+&
                    (1d0-theta)*matmul(transpose(Sk_nc_old),rk_old))
            end if
            if (mix_waters(j+1)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
                Se_nc=mix_waters(j+1)%solid_chemistry%reactive_zone%get_Se_nc_react_zone(ind_var_act_species_rz) !> we get the equilibrium stoichiometric matrix in reactive zone
                allocate(U_SeT(size(comp_mat_rz,1),size(Se_nc,1)))
                U_SeT=matmul(comp_mat_rz,transpose(Se_nc)) !> we compute U*Se_nc^T
                if (norm_mat_inf(U_SeT)>=tol) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Se_nc),&
                        mix_waters(j+1)%solid_chemistry%re_mean)
                end if
                deallocate(U_SeT)
            end if
        end do
        Sk_nc_up=Sk_nc !> chapuza
        Sk_nc_up_old=Sk_nc_old !> chapuza
        rk_up_old=rk_old !> chapuza
    !> Downstream waters
        do j=1,num_down
            if (mix_waters(num_up+1+j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+1+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                Sk_nc_old=mix_waters(num_up+1+j)%get_Sk_nc_old() !> chapuza
                if (norm_mat_inf(Sk_nc_up-Sk_nc)<tol) then !> we check if the kinetic stoichiometric matrix is the same
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+1+j)*(&
                        theta*matmul(transpose(Sk_nc_up),rk_new)+&
                        (1d0-theta)*matmul(transpose(Sk_nc_up_old),rk_up_old))
                else
                    rk_old=mix_waters(num_up+1+j)%get_rk_old() !> chapuza
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+1+j)*&
                        matmul(transpose(Sk_nc_old),rk_old)
                end if
            end if
            if (mix_waters(num_up+1+j)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
                Se_nc=mix_waters(num_up+1+j)%solid_chemistry%reactive_zone%get_Se_nc_react_zone(ind_var_act_species_rz) !> we get the equilibrium stoichiometric matrix in reactive zone
                allocate(U_SeT(size(comp_mat_rz,1),size(Se_nc,1)))
                U_SeT=matmul(comp_mat_rz,transpose(Se_nc)) !> we compute U*Se_nc^T
                if (norm_mat_inf(U_SeT)>=tol) then
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+1+j)*matmul(transpose(Se_nc),&
                        mix_waters(num_up+1+j)%get_re_mean())
                end if
                deallocate(U_SeT)
            end if
        end do
    !> Current water
    if (mix_waters(1)%indices_rk%num_cols>0) then
        Sk_nc_old=mix_waters(1)%get_Sk_nc_old() !> we get the kinetic stoichiometric matrix in 1-th upstream water (chapuza)
        rk_old=mix_waters(1)%get_rk_old() !> chapuza
        r_tilde=r_tilde+mixing_ratios_reacts(1)*(&
            (1d0-theta)*matmul(transpose(Sk_nc_old),rk_old))
    end if
end subroutine

!> This subroutine computes the kinetic reaction rate contributions after water mixing implicitly for both upstream and downstream waters
!> The estimation of the downstream waters reaction amoutns is done by extrapolation with the second time derivative of the reaction rates
!> We assume uniform time step
subroutine compute_r_tilde_impl_opt3(mix_waters,ind_var_act_species_rz,comp_mat_rz,&
    mixing_ratios_reacts,num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other (must be ordered in current, upstream & downstream waters, respectively)
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> indices of variable activity species in the current reactive zone
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
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: conc_nc(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: rk_est(:) !> auxiliary variable to store the kinetic reaction rates after mixing
    !logical :: flag !> flag to check if an index is in "perm"
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    !real(kind=8) :: mu=5d-1 !> parameter for extrapolation of downstream waters rk
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    
    !allocate(c_mix(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    !allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down+1) error stop "Dimension error in compute_r_tilde_impl_opt3"
    
!> Process
    !c_mix=0d0 !> initialise concentration of aqueous variable activity species after mixing
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    r_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j+1)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_R(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Sk_nc),&
                        theta*mix_waters(j+1)%rk_new+&
                        (1d0-theta)*mix_waters(j+1)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j+1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    r_tilde=r_tilde+mixing_ratios_reacts(j+1)*matmul(transpose(Sk_nc),&
                        theta*mix_waters(j+1)%solid_chemistry%rk_new+&
                        (1d0-theta)*mix_waters(j+1)%solid_chemistry%rk_old)
                end if
                deallocate(Sk_nc)
            end if
        end do
    !> Downstream waters (we estimate them)
        do j=1,num_down
            if (mix_waters(num_up+j+1)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+j+1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                if (mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    allocate(rk_est(mix_waters(num_up+j+1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts))
                    rk_est=(1d0+&
                        2d0*theta)*mix_waters(num_up+j+1)%rk_old-&
                        3d0*theta*mix_waters(num_up+j+1)%rk_old_old+&
                        theta*mix_waters(num_up+j+1)%rk_old_old_old
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j+1)*matmul(transpose(Sk_nc),&
                        rk_est)
                    call mix_waters(num_up+j+1)%set_Rk_est(rk_est*Delta_t)
                    !end if
                    deallocate(rk_est)
                end if
                if (allocated(mix_waters(num_up+j+1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    allocate(rk_est(mix_waters(num_up+j+1)%solid_chemistry%mineral_zone%num_minerals_kin))
                    rk_est=(1d0+&
                        2d0*theta)*mix_waters(num_up+j+1)%solid_chemistry%rk_old-&
                        3d0*theta*mix_waters(num_up+j+1)%solid_chemistry%rk_old_old+&
                        theta*mix_waters(num_up+j+1)%solid_chemistry%rk_old_old_old
                    r_tilde=r_tilde+mixing_ratios_reacts(num_up+j+1)*matmul(transpose(Sk_nc),&
                        rk_est)
                    !end if
                    call mix_waters(num_up+j+1)%solid_chemistry%set_Rk_est(rk_est*Delta_t)
                    deallocate(rk_est)
                end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
        !> Current water
        if (mix_waters(1)%indices_rk%num_cols>0) then
            Sk_nc=mix_waters(1)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in 1-th upstream water (chapuza)
            !rk_tilde=rk_tilde+mixing_ratios_R(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
            !deallocate(Sk_nc)
            if (mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc),&
                    mix_waters(1)%rk_old)
                !deallocate(Sk_nc)
            end if
            if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                r_tilde=r_tilde+(1d0-theta)*mixing_ratios_reacts(1)*matmul(transpose(Sk_nc),&
                    mix_waters(1)%solid_chemistry%rk_old)
            end if
            !Sk_nc_old=Sk_nc !> chapuza
            deallocate(Sk_nc)
        end if
end subroutine 

!> This subroutine computes the kinetic reaction rate contributions after water mixing, explicitly for all target waters
subroutine compute_r_tilde_expl(mix_waters,ind_var_act_species_rz,comp_mat_rz,mixing_ratios_reacts,&
    num_up,num_down,theta,Delta_t,r_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    integer(kind=4), intent(in) :: ind_var_act_species_rz(:)  !> indices of variable activity species in the current reactive zone
    real(kind=8), intent(in) :: comp_mat_rz(:,:)  !> indices of variable activity species in the current reactive zone
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
                r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk_old)
            end if
            if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                r_tilde=r_tilde+mixing_ratios_reacts(j)*matmul(transpose(Sk_nc),mix_waters(j)%solid_chemistry%rk_old)
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
!! \details This subroutine establishes the index mappings between aqueous species in the
!! aqueous_chemistry object and the corresponding species in the chemical system. It creates
!! three key index arrays:
!! - **ind_prim_species**: Indices of primary (basis) aqueous species
!! - **ind_sec_species**: Indices of secondary aqueous species (complexes + constant activity species)
!! - **ind_var_act_species**: Indices of all variable activity species (primary + secondary variable)
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
            
            !> \subsection alt_implementation Alternative implementation (commented out)
            !! This commented code represents an alternative approach that was considered but not used.
            !! It would iterate through all species and classify them based on their position
            !! relative to num_aq_prim_species, but the current implementation above is preferred
            !! because it explicitly iterates through primary and secondary species separately.
            
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
            
            !> \subsection alt_aux_array Alternative auxiliary array (commented out)
            !! This commented code would use an auxiliary array to temporarily store indices
            !! before rearranging them. Not currently used.
            
            !> Deallocate existing arrays if allocated (ensures clean state)
            
            !> Allocate index arrays based on speciation algebra dimensions
            !allocate(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
            
            !> \name Loop: Process all species with activity-based classification
            !! Iterate through all aqueous species and classify them based on three criteria:
            !! 1. Constant activity → goes to end of ind_sec_species
            !! 2. Variable activity in "primary" positions → goes to ind_prim_species and ind_var_act_species
            !! 3. Variable activity in "secondary" positions → goes to ind_sec_species and ind_var_act_species
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
        !! OJO: FALTAN REACCIONES CINÉTICAS LINEALES
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
        rk_new(i)=this%rk_new(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rk_new(i)
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
        r_new(n_r+i)=this%solid_chemistry%re_mean(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%num_minerals_cst_act
    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        r_new(n_r+i)=this%gas_chemistry%re_mean(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
        r_new(n_r+i)=this%re_mean(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        r_new(n_r+i)=this%rk_new(i)
    end do
    n_r=n_r+this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        r_new(n_r+i)=this%solid_chemistry%rk_new(i)
    end do
    !> Faltan gases
    end function

    function get_r_old(this) result(r_old)
    !!> This function returns the old kinetic reaction rates ordered in aqueous & mineral reactions (ESTA MAL)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: r_old(:)
    integer(kind=4) :: i, n_r
    allocate(r_old(this%indices_rk_old%num_cols+this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions))
    n_r=0
    do i=1,this%solid_chemistry_old%reactive_zone%num_minerals_cst_act
        r_old(n_r+i)=this%solid_chemistry_old%re_mean(i)
    end do
    n_r=n_r+this%solid_chemistry_old%reactive_zone%num_minerals_cst_act
    do i=1,this%solid_chemistry_old%reactive_zone%gas_phase%num_gases_eq_cst_act
        r_old(n_r+i)=this%gas_chemistry%re_mean(i) !> tiene que ser de gas_chemistry_old
    end do
    n_r=n_r+this%solid_chemistry_old%reactive_zone%gas_phase%num_gases_eq_cst_act
    do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts
        r_old(n_r+i)=this%re_mean(i)
    end do
    n_r=n_r+this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_eq_reacts
    do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts
        r_old(n_r+i)=this%rk_old(i)
    end do
    n_r=n_r+this%solid_chemistry_old%reactive_zone%chem_syst%num_aq_kin_reacts
    do i=1,this%solid_chemistry_old%mineral_zone%num_minerals_kin
        r_old(n_r+i)=this%solid_chemistry_old%rk_old(i)
    end do
    !> Faltan gases
    end function
    
    function get_Rk(this) result(Rk)
    !!> This function returns the kinetic reaction amounts ordered in aqueous & mineral reactions
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: Rk(:)
    integer(kind=4) :: i
    allocate(Rk(this%indices_rk%num_cols))
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        Rk(i)=this%Rk(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        Rk(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%Rk(i)
    end do
    !> Faltan gases
    end function

    function get_rk_old(this) result(rk_old)
        !!> This function returns the old kinetic reaction rates (chapuza) ordered in aqueous & mineral reactions
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: rk_old(:)
        integer(kind=4) :: i
        allocate(rk_old(this%indices_rk%num_cols))
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
            rk_old(i)=this%rk_old(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rk_old(i)
        end do
        !> Faltan gases
    end function

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
        integer(kind=4), intent(in), optional :: ind_var_act_species(:) !< Indices of variable activity species in chemical system
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

    !> \brief Returns the variable activity stoichiometric matrix from the previous time step
    !!
    !! @param ind_var_act_species optional indices of variable activity species
    !! @return S_nc_old stoichiometric matrix from previous time step
    !!
    function get_S_nc_old(this,ind_var_act_species) result(S_nc_old) !> returns variable activity stoichiometric matrix of kinetic reactions
        class(aqueous_chemistry_c) :: this
        integer(kind=4), intent(in), optional :: ind_var_act_species(:) !< Indices of variable activity species in chemical system
        real(kind=8), allocatable :: S_nc_old(:,:)
        integer(kind=4) :: i,n_r,num_rk, j
        allocate(S_nc_old(this%indices_rk_old%num_cols+this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions,&
            this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_Act_species))
        ! S_nc_old(1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions,:)=&
        !     this%solid_chemistry_old%reactive_zone%get_Se_nc_react_zone(ind_var_act_species)
        do i=1,this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions     
            do j=1,size(ind_var_act_species)
                !< Extract column for each specified species
                S_nc_old(i,j)=this%solid_chemistry_old%reactive_zone%stoich_mat(&
                    i,this%solid_chemistry_old%reactive_zone%ind_var_act_species(ind_var_act_species(j)))
            end do
        end do
        n_r=this%solid_chemistry_old%reactive_zone%speciation_alg%num_eq_reactions
        num_rk=0
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
            S_nc_old(n_r+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
        n_r=n_r+num_rk
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
            S_nc_old(n_r+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        num_rk=num_rk+this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
        n_r=n_r+num_rk
        do i=1,this%solid_chemistry_old%mineral_zone%num_minerals_kin
            S_nc_old(n_r+i,:)=&
                this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+&
                this%solid_chemistry_old%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        !> Faltan gases
    end function

    !> \brief Returns the kinetic stoichiometric matrix from the previous time step
    !!
    !! @return Sk_nc_old kinetic stoichiometric matrix from previous time step
    !!
    function get_Sk_nc_old(this) result(Sk_nc_old) !> returns variable activity stoichiometric matrix of kinetic reactions (chapuza)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: Sk_nc_old(:,:)
        integer(kind=4) :: i,num_rk
        allocate(Sk_nc_old(this%indices_rk_old%num_cols,this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species))
        num_rk=0
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
            Sk_nc_old(num_rk+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        end do
        num_rk=num_rk+this%solid_chemistry_old%reactive_zone%chem_syst%num_lin_kin_reacts
        do i=1,this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
            Sk_nc_old(num_rk+i,:)=this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+i,&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        end do
        num_rk=num_rk+this%solid_chemistry_old%reactive_zone%chem_syst%num_redox_kin_reacts
        do i=1,this%solid_chemistry_old%mineral_zone%num_minerals_kin
            Sk_nc_old(num_rk+i,:)=&
                this%solid_chemistry_old%reactive_zone%chem_syst%Sk(num_rk+&
                this%solid_chemistry_old%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry_old%reactive_zone%speciation_alg%num_var_act_species)
        end do
        
        !> Faltan gases
    end function
    
    !> \brief Computes kinetic reaction amounts during a time step
    !!
    !! \details Uses time-weighted average: Rk = (theta*rk_new + (1-theta)*rk_old)*Delta_t
    !!
    !! @param theta time weighting factor
    !! @param Delta_t time step
    !!
    subroutine compute_Rk(this,theta,Delta_t) !> computes kinetic reaction amounts during a time step
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    integer(kind=4) :: i
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        this%Rk(i)=theta*this%rk_new(i)+(1d0-theta)*this%rk_old(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        this%solid_chemistry%Rk(i)=theta*this%solid_chemistry%rk_new(i)+(1d0-theta)*this%solid_chemistry%rk_old(i)
    end do
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        this%Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
            theta*this%rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)+&
            (1d0-theta)*this%rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
    end do
    this%Rk=this%Rk*Delta_t
    this%solid_chemistry%Rk=this%solid_chemistry%Rk*Delta_t
    end subroutine

    !> \brief Sets the mean kinetic reaction rates during a time step
    !!
    !! @param rk_mean mean kinetic reaction rates during a time step
    !!
    subroutine set_rk_mean(this,rk_mean)
        !!> This subroutine sets the mean kinetic reaction rates during a time step
        !! We assume that rk_mean is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: rk_mean(:) !> mean kinetic reaction rates during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%rk_mean(i)=rk_mean(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rk_mean(i)=rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
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
            this%Rk(i)=Rk(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%Rk(i)=Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
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
            this%rk_old(i)=rk_old(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rk_old(i)=rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)
        end do
        !> Faltan gases
    end subroutine
    
    !> \brief Sets the estimated kinetic reaction rates
    !!
    !! \details Sets the estimated kinetic reaction rates in the consistent
    !! WMA approach.
    !!
    !! @param Rk_est estimated kinetic reaction rates
    !!
    subroutine set_Rk_est(this,Rk_est)
    !!> This subroutine sets the estimated kinetic reaction rates in the consistent WMA approach
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: Rk_est(:) !> estimated kinetic reaction rates
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%Rk_est(i)=Rk_est(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%Rk_est(i)=Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
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
        if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
            flag_surf=.true.
        else
            flag_surf=.false.
            !> chapuza: to ensure that solid chemistry arrays are allocated
            call this%solid_chemistry%allocate_conc_solids()
            call this%solid_chemistry%allocate_activities()
            call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
        end if
        call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
        call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions()
        call this%set_ind_species()
        call this%solid_chemistry%reactive_zone%set_ind_eq_reacts(ind_eq_reacts)
        call this%solid_chemistry%reactive_zone%allocate_ind_var_act_species()
        call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
        call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
        call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
        !allocate(swap(2))
        if (associated(this%gas_chemistry)) then
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap,this%gas_chemistry%activities)
        else
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
        end if
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
        if (associated(this%gas_chemistry)) then
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap,this%gas_chemistry%activities)
        else
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
        end if
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
        this%re_mean=this%Re/Delta_t
        this%solid_chemistry%re_mean=this%solid_chemistry%Re/Delta_t
        this%gas_chemistry%re_mean=this%gas_chemistry%Re/Delta_t
    end subroutine

    !> Compute the concentration of ideal water in the solution
    !! Mathematical Context: Water molality = 1000 g / 18.015 g/mol ≈ 55.51 mol/kg H2O
    !! For ideal dilute solutions, water activity ≈ 1, and this concentration is used for activity calculations
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
        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        
        !> Copy scalar position array if allocated in source
        
        !> Copy all scalar thermodynamic properties
        
        !> Set pointer association to solid chemistry object (not deep copy)
        if (associated(other%solid_chemistry)) then                        !< Check if solid chemistry pointer is associated in source
            call this%set_solid_chemistry(other%solid_chemistry)           !< Set pointer association to same solid chemistry object
        end if
        
        !> Set pointer association to gas chemistry object (not deep copy)
        if (associated(other%gas_chemistry)) then                          !< Check if gas chemistry pointer is associated in source
            call this%set_gas_chemistry(other%gas_chemistry)               !< Set pointer association to same gas chemistry object
        end if
        
        !> Set pointer association to aqueous phase object (not deep copy)
        
        !> Deep copy all index arrays (independent copies needed for different chemistry objects)
        if (allocated(other%ind_prim_species)) then                        !< Check if primary species indices are allocated in source
            this%ind_prim_species=other%ind_prim_species                   !< Deep copy primary species indices array
        end if
         if (allocated(other%ind_var_act_species)) then                     !< Check if variable activity species indices are allocated in source
             this%ind_var_act_species=other%ind_var_act_species             !< Deep copy variable activity species indices array
         end if
        if (allocated(other%ind_sec_species)) then                         !< Check if secondary species indices are allocated in source
            this%ind_sec_species=other%ind_sec_species                     !< Deep copy secondary species indices array
        end if
        this%indices_rk=other%indices_rk                                   !< Copy reaction kinetics indices structure (value assignment)
    end subroutine                                                         !< End of assignment subroutine

    subroutine copy_aq_chem_Lagr(this,other)
        !!> This subroutine assigns the content of an aqueous chemistry object to another one
        class(aqueous_chemistry_c) :: this                                 !< Destination aqueous chemistry object (target for assignment)
        class(aqueous_chemistry_c), intent(in) :: other                    !< Source aqueous chemistry object (source for copying)

        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        call this%copy_local_chemistry(other)                            !< Assign parent class data from local_chemistry_c base class
        !print *, associated(other%gas_chemistry)                          !< Debug print statement (commented out)
        
        !> Copy scalar position array if allocated in source
        !if (allocated(other%pos)) then                                     !< Check if position array is allocated in source object
        
        !> Copy all scalar thermodynamic properties
        this%params_aq_sol=other%params_aq_sol                             !< Copy aqueous solution parameters structure
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
    end subroutine    
    
        !> \brief Sets indices of aqueous variable activity species in the state variables
        subroutine set_ind_aq_var_act_species(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object

            integer(kind=4) :: i !> loop indices

            if (allocated(this%ind_aq_var_act_species)) then
                deallocate(this%ind_aq_var_act_species) !> deallocate if already allocated
            end if

            allocate(this%ind_aq_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !> allocate array for target waters

            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> for each aqueous primary species
                this%ind_aq_var_act_species(i)=this%ind_var_act_species(i) !> set index
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !> for each aqueous secondary variable activity species
                this%ind_aq_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)=&
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i) !> set index
            end do
        end subroutine
        
        !> \brief Computes logarithmic activity coefficients for the aqueous phase
        !!
        !! \details Calculates ionic strength and then uses aqueous solution
        !! parameters to compute log activity coefficients for all aqueous species.
        !!
        subroutine compute_log_act_coeffs_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            
            call this%compute_ionic_strength()
            call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_strength,this%params_aq_sol,this%log_act_coeffs)
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
            d_log_gamma_d_I=this%compute_d_log_gamma_d_I_aq_chem() !> Compute ∂(log γᵢ)/∂I for all variable activity species (from Debye-Hückel or extended models)
            out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:n_v_aq)) !> Compute outer product: [∂(log γ)/∂I] ⊗ [z²] for activity corrections in Jacobian
            call this%aq_phase%compute_log_Jacobian_act_coeffs_aq_phase(out_prod,this%concentrations(this%ind_aq_var_act_species),this%log_Jacobian_act_coeffs(1:n_v_aq,1:n_v_aq))
            deallocate(d_log_gamma_d_I)
        end subroutine

        !> \brief Get logarithmic activity coefficients of all secondary species.
        !>
        !> This function retrieves the logarithmic activity coefficients (log γ₂) for all secondary species
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
        !> \return log_gamma2 Vector of log γ₂ values [-] (size = num_sec_species)
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
            !!> Kinetic variable activity minerals
            !!> Kinetic variable activity gases
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
            ! !> Ideal water
            ! !> Kinetic constant activity minerals
            
            ! !> Kinetic constant activity gases
            !> Equilibrium constant activity gases
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                log_gamma2(n2+i)=this%gas_chemistry%log_act_coeffs(i)
            end do
            n2=n2+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
        end function

        !> \brief Get logarithmic activity coefficients of all primary species.
        !>
        !> This function retrieves the logarithmic activity coefficients (log₁₀ γ₁) for all primary
        !> species in the chemical system. Primary species (also called basis species or components)
        !> are the fundamental independent species used as a basis for expressing all other species
        !> through mass action relationships in chemical equilibrium calculations.
        !>
        !> The log activity coefficients are assembled in the following order:
        !> 1. Aqueous primary species (basis species for aqueous complexation, e.g., Ca²⁺, HCO₃⁻)
        !> 2. Free sites (unoccupied surface sites for cation exchange, if present)
        !> 3. Kinetic variable activity minerals (minerals undergoing kinetic dissolution/precipitation)
        !> 4. Kinetic variable activity gases (gases with time-dependent partial pressures)
        !>
        !> \param[in] this Aqueous chemistry object with current chemical state and activity coefficients
        !> \param[out] log_gamma1 Logarithmic activity coefficients of all primary species (allocated internally)
        !> \return log_gamma1 Vector of log₁₀ γ₁ values [-] (size = num_prim_species)
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
                log_gamma_var_act(i)=0d0 !> ideal solution: log γ = 0
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
                log_gamma_var_act(n_v+i)=0d0 !> ideal solution: log γ = 0
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
        
!> Modify mixing ratios for kinetic reaction rates to prevent negative concentrations
!! Mathematical Context: Adaptive time stepping for reactive transport to maintain positivity
!! Algorithm: If c + Δt×r < 0, reduce mixing ratio λ to ensure c + Δt×λ×r ≥ 0
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

        rk_old=this%get_rk_old()                                   !< Get old kinetic reaction rates
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
                mix_ratios_new(1)=mix_ratios_new(1)*alpha                    !< Increase mixing ratio by factor α
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
                ! mix_ratios_new(1)=mix_ratios_init(1)*alpha                    !< Increase mixing ratio by factor α
                !mix_ratios_new(2:)=mix_ratios_new(2:)*beta                    !< Increase mixing ratio by factor α
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
        rk_old=this%get_rk_old()                                   !< Get old kinetic reaction rates
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
                re_mean(n_eq+i)=this%solid_chemistry%re_mean(i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%num_minerals_cst_act
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                re_mean(n_eq+i)=this%gas_chemistry%re_mean(i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
                re_mean(n_eq+i)=this%re_mean(n_eq+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
            do i=1,this%solid_chemistry%reactive_zone%num_minerals_var_act
                re_mean(n_eq+i)=this%solid_chemistry%re_mean(&
                    this%solid_chemistry%reactive_zone%num_minerals_cst_act+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%num_minerals_var_act
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                re_mean(n_eq+i)=this%solid_chemistry%re_mean(&
                    this%solid_chemistry%reactive_zone%num_minerals+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                re_mean(n_eq+i)=this%gas_chemistry%re_mean(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
            end do
            n_eq=n_eq+this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
        end function get_re_mean

        subroutine update_old_attributes(this)
        !> This subroutine updates the old attributes to the current state
            class(aqueous_chemistry_c), intent(inout) :: this                       !< Aqueous chemistry object (input/output)
            call this%update_conc_old() !< Update old log activity coefficients
            call this%solid_chemistry%update_conc_old() !< Update old solid chemistry state
            !call this%gas_chemistry%update_conc_old() !< Update old gas chemistry state
            call this%update_rk_old() !< Update old reaction rates
            call this%solid_chemistry%update_rk_old() !< Update indices for old reaction rates
            this%indices_rk_old=this%indices_rk !< Update indices for old reaction rates
            this%solid_chemistry_old=>this%solid_chemistry !< Update old solid chemistry state
            this%ind_prim_species_old=this%ind_prim_species !< Update old primary species indices
            this%ind_sec_species_old=this%ind_sec_species !< Update old secondary species indices
        end subroutine update_old_attributes
end module                                                                  !< End of aqueous_chemistry_m module
