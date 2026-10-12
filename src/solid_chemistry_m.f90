!> \file solid_chemistry_m.f90
!> \brief Module for solid-phase chemistry, including minerals, surface complexation, and cation exchange.
!>
!> This module defines the `solid_chemistry_c` type and its associated procedures for managing solid-phase chemical processes in reactive transport simulations. It supports mineral precipitation/dissolution, surface complexation, cation exchange, and coupling to aqueous chemistry. The module provides memory management, property calculation, iterative solvers, and validation routines for solid-phase chemical state.
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2025
!> \version 1.0
module solid_chemistry_m
    use local_chemistry_m, only: local_chemistry_c                         !< Import base chemistry class for inheritance
    use reactive_zone_m, only: reactive_zone_c
    use vectors_m, only: inf_norm_vec_real
    use biofilm_m, only: biofilm_c                                 !< Import biofilm class for biofilm properties
    use mineral_zone_m, only: mineral_zone_c                    !< Import mineral zone and mineral classes
    use target_m, only: target_c                                           !< Import target class for spatial discretization
    implicit none                                                           !< No implicit variable declarations
    save                                                                    !< Preserve module variables between procedure calls
    private                                                                 !< Default accessibility is private
!***********************************************************************************************************************************
    !> Solid chemistry class: Manages all solid phase chemical processes
    !! Inheritance: Extends local_chemistry_c to add solid-specific functionality
    !! Design Pattern: Composition pattern with pointers to reactive and mineral zones
    type, public, extends(local_chemistry_c) :: solid_chemistry_c
        !=====================================================================
        ! COMMON / SHARED SOLID-PHASE STATE
        !=====================================================================
        integer(kind=4) :: num_solids=0                                     !< Total number of solid species (temporary variable - "chapucilla")
        class(reactive_zone_c), pointer :: reactive_zone                   !< Pointer to reactive zone containing equilibrium solids (minerals, surface complexes, gases)
        class(target_c), pointer :: tar                                    !< Pointer to associated spatial target (bijective mapping)
        real(kind=8) :: time                       !< Species travel time through solid portion, used for kinetic rate calculations (size = num_solids)
        logical :: is_interface = .false.                           !< Flag indicating whether the solid chemistry object is an interface
        !=====================================================================
        ! MINERALS
        !=====================================================================
        class(mineral_zone_c), pointer :: mineral_zone                     !< Pointer to mineral zone containing kinetic and equilibrium minerals
        real(kind=8), allocatable :: min_pore_vol_ratios(:)                         !< Mineral-to-pore volume ratios V_min/V_pore (dimensionless, >= 0, may exceed 1)
                                                                            !! Order: kinetic minerals first, then equilibrium minerals
        real(kind=8), allocatable :: react_surfaces(:)                     !< Reactive surface area of each mineral per unit volume of water [m² mineral/L water] as read from the input (aqueous_chemistry_c%compute_react_surfaces_per_kg_water returns the per-kg-H2O values without modifying this array)
                                                                            !! Used for surface complexation and mineral kinetics

        !=====================================================================
        ! CATION EXCHANGE
        !=====================================================================
        real(kind=8) :: CEC                                                !< Cation Exchange Capacity [eq/L_bulk]
                                                                            !! Assumption: maximum one adsorption surface per solid chemistry object
        real(kind=8), allocatable :: equivalents(:)                        !< Equivalents of adsorbed cations [eq/L_bulk] for charge balance

        !=====================================================================
        ! BIOFILM
        !=====================================================================
        class(biofilm_c), pointer :: biofilm                               !< Pointer to associated biofilm object (intrinsic, medium-independent parameters)
        real(kind=8), allocatable :: vol_fracts_microorgs(:)                !< Microorganism volumetric fractions: vol_fracs(i) = V_microorg_i / V_total_biomass (size = num_microorgs) [-]
        !> Biofilm LOCAL STATE (medium-dependent; the intrinsic biofilm parameters live in the biofilm pointer above)
        real(kind=8) :: f_w                          !< Water volumetric fraction in biofilm: f_w = V_water / V_biofilm [-]
        real(kind=8) :: f_cells                      !< Cell volume fraction in biofilm [-]
        real(kind=8) :: f_eps                        !< EPS volume fraction in biofilm [-]
        real(kind=8) :: thickness                    !< Biofilm thickness Lf [m]
        real(kind=8) :: o2_penetration_depth         !< O2 penetration depth [m]
        real(kind=8) :: bl_thickness                 !< External mass-transfer boundary layer thickness L_bl [m]
        real(kind=8) :: roughness_Ra                 !< Surface roughness coefficient Ra [-]
        real(kind=8) :: fractal_dim_2d               !< 2D fractal dimension Df [-]
        real(kind=8) :: fractal_dim_3d               !< 3D fractal dimension Df [-]
        real(kind=8) :: surface_coverage             !< Surface coverage fraction [-]
        real(kind=8) :: specific_biovolume           !< Specific biovolume per substratum area [m^3/m^2]
        real(kind=8) :: porosity_initial             !< Initial (clean) medium porosity n_0 [-]
        real(kind=8) :: porosity_biofilm             !< Biofilm-affected medium porosity n_b [-]
        real(kind=8) :: pH_min                       !< Minimum interior pH [-]
        real(kind=8) :: pH_max                       !< Maximum interior pH [-]
        real(kind=8) :: Eh_aerobic                   !< Redox potential in aerobic zone [mV vs SHE]
        real(kind=8) :: Eh_anaerobic                 !< Redox potential in anaerobic core [mV vs SHE]
    contains
        procedure :: initialise_solid_chemistry                    !< Compute the number of solids and allocate every state array from the associated reactive and mineral zones
    !> Setter procedures - Assign values to solid chemistry properties
        procedure :: set_concentrations=>set_conc_solids           !< Set solid species concentrations [mol/L_bulk]
        procedure :: set_indices_solids                            !< Set index mapping for solid species ordering
        procedure :: set_vol_fracts                                !< Set volumetric fractions of minerals
        procedure :: set_react_surfaces                            !< Set specific surface areas for minerals
        procedure :: set_reactive_zone                             !< Associate with reactive zone object
        procedure :: set_mineral_zone                              !< Associate with mineral zone object
        procedure :: set_CEC                                       !< Set cation exchange capacity
        procedure :: set_conc_free_site                            !< Set concentration of free surface sites
        procedure :: set_act_surf_compl                            !< Set activities of surface complexes
        procedure :: set_target                                    !< Associate with spatial target object
        procedure :: set_time                                    !< Associate with spatial target object
        procedure :: set_is_interface                              !< Set the interface flag for the solid chemistry object
    !> Memory allocation procedures - Allocate arrays for solid chemistry data
        procedure :: allocate_min_pore_vol_ratios                     !< Allocate volumetric fractions array
        procedure :: allocate_vol_fracts_microorgs                !< Allocate volumetric fractions array
        procedure :: allocate_react_surfaces                       !< Allocate reactive surfaces array
        procedure :: allocate_conc_solids                          !< Allocate solid concentrations array
        procedure :: allocate_rxn_rates_and_amounts_solid_chem            !< Allocate reaction rates array for solid reactions
        procedure :: deallocate_rxn_rates_and_amounts_solid_chem          !< Release the reaction rate/amount arrays (re-layout for new zones)
        procedure :: allocate_activities                           !< Allocate activities array for solid species
        procedure :: allocate_log_act_coeffs_solid_chem            !< Allocate log activity coefficients array
        procedure :: allocate_equivalents                          !< Allocate equivalents array for charge balance
    !> Computational procedures - Calculate solid chemistry properties and equilibrium
        procedure :: compute_activities_solids                     !< Compute activities of solid species (usually = concentrations)
        procedure :: compute_conc_surf_compl                       !< Compute surface complex concentrations from activities and activity coefficients
        procedure :: compute_equivalents                           !< Compute equivalents of adsorbed cations for charge balance
        procedure :: compute_conc_minerals_iter                    !< Iteratively compute mineral concentrations for equilibrium
        procedure :: compute_mass_bal_mins                         !< Compute mineral mass balance errors for validation
        procedure :: compute_conc_surf_ideal_Newton                !< Compute surface complexation using Newton-Raphson method
        procedure :: compute_conc_surf_ideal_Picard                !< Compute surface complexation using Picard iteration
        procedure :: compute_conc_surf_ideal_bin                   !< Compute surface complexation assuming ideal behavior for binary cation exchange systems
        procedure :: compute_num_solids_solid_chem                 !< Count total number of solid species
    !> Biofilm procedures - Compute biofilm properties from local state + biofilm parameters
        procedure :: compute_permeability_reduction                !< Kozeny-Carman bioclogging permeability ratio [-]
        procedure :: compute_detachment_rate                       !< Biofilm detachment rate [m/s]
        procedure :: compute_thiele_modulus                        !< Biofilm Thiele modulus phi^2 [-]
    !> Getter procedures - Retrieve solid chemistry properties
        procedure :: get_tar_id                                     !< Get target ID for solid chemistry object
    !> Update procedures - Modify solid chemistry state during simulations
        procedure :: update_conc_ads_cats                          !< Update concentrations of adsorbed cations
        procedure :: update_act_ads_cats                           !< Update activities of adsorbed cations
    !> Validation procedures - Check consistency and validity of solid chemistry state
        procedure :: check_solid_chemistry                         !< Validate solid chemistry object for consistency
    !> Assignment procedures - Copy and manipulate solid chemistry objects
        procedure :: copy_solid_chemistry                        !< Deep copy assignment of solid chemistry object
        procedure :: modify_mix_ratios_reacts_solid_chem                          !< Modify mixing ratios for reaction amounts
    end type
    
    !> Interface block for surface complexation Newton-Raphson solver
    !! Mathematical Context: Newton-Raphson method for nonlinear surface complexation equilibrium
    !! Solves: Σ(Xᵢ) = CEC where Xᵢ are adsorbed cation concentrations
    interface
        subroutine compute_conc_surf_ideal_Newton(this,conc_cats,act_ads_cats_ig,niter,CV_flag,logK)
            import solid_chemistry_c                                        !< Import solid chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_cats(:)                       !< Cation concentrations in solution (temporary - "chapuza")
                                                                            !! Dimension: number of cation exchange half reactions
            real(kind=8), intent(in) :: act_ads_cats_ig(:)                 !< Initial guess for adsorbed cation activities (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
            real(kind=8), intent(in) :: logK(:)                            !< Per-water log10 equilibrium constants of the EXCHANGE reactions only (size = num_exch_cats), slice of aqueous_chemistry_c%logK
        end subroutine
        
        !> Interface for surface complexation Picard iteration solver
        !! Mathematical Context: Picard iteration (fixed-point) method for surface complexation
        !! More robust than Newton-Raphson but slower convergence rate
        subroutine compute_conc_surf_ideal_Picard(this,conc_cats,act_ads_cats_ig,niter,CV_flag)
            import solid_chemistry_c                                        !< Import solid chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_cats(:)                       !< Cation concentrations in solution (temporary - "chapuza")
                                                                            !! Dimension: number of cation exchange half reactions
            real(kind=8), intent(in) :: act_ads_cats_ig(:)                 !< Initial guess for adsorbed cation activities (pre-allocated)
            integer(kind=4), intent(out) :: niter                          !< Number of iterations performed
            logical, intent(out) :: CV_flag                                !< Convergence flag: TRUE if converges, FALSE otherwise
        end subroutine
        
        !> Interface for direct surface complexation calculation (ideal case)
        !! Mathematical Context: Direct solution assuming ideal behavior (no activity corrections)
        !! Used when surface complexation is weak or as initial approximation
        subroutine compute_conc_surf_ideal_bin(this,conc_cats,logK)
            import solid_chemistry_c                                        !< Import solid chemistry class
            implicit none                                                   !< No implicit variable declarations
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc_cats(:)                       !< Cation concentrations in solution (temporary - "chapuza")
                                                                            !! Dimension: number of cation exchange half reactions
            real(kind=8), intent(in) :: logK(:)                            !< Per-water log10 equilibrium constants of the EXCHANGE reactions only (size = num_exch_cats), slice of aqueous_chemistry_c%logK
            !integer(kind=4), intent(out) :: niter                         !< Commented iteration counter
            !logical, intent(out) :: CV_flag                               !< Commented convergence flag
        end subroutine
        
      
        
        subroutine compute_mass_bal_mins(this,Delta_t)
            import solid_chemistry_c
            implicit none
            class(solid_chemistry_c) :: this
            !real(kind=8), intent(in) :: re_mean(:) !> equilibrium reaction rates
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio
        end subroutine
        
       
    end interface                                                               !< End of interface block
    
    
    
    contains                                                                !< Beginning of implementation section
        
       
        
        !> Set concentrations of solid species
        !! Physical Context: Solid concentrations typically in [mol/L_bulk] including both minerals and surface complexes
        !! Validation: Ensures input array size matches number of solid species in chemical system
        subroutine set_conc_solids(this,conc)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: conc(:)                            !< Solid species concentrations array [mol/L_bulk]
            if (this%reactive_zone%chem_syst%num_solids<size(conc)) error stop "Dimension error in set_conc_solids"
                                                                            !< Validate array size against number of solid species
            this%concentrations=conc                                       !< Assign concentrations to object
        end subroutine        
        
        !> Set concentration of free surface sites to default value
        !! Physical Context: Free sites are available for surface complexation reactions
        !! Default value: 1e-32 mol/L_bulk (very small but non-zero to avoid numerical issues)
        !! Logic: Handles both mineral zone and reactive zone associations
        subroutine set_conc_free_site(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            !real(kind=8), intent(in) :: conc                              !< Commented concentration parameter
            if (ASSOCIATED(this%mineral_zone)) then                        !< If mineral zone is associated
                this%concentrations(this%mineral_zone%num_minerals+1)=1d-32   !< Set free site concentration after minerals
                                                                            !! Position: after kinetic and equilibrium minerals
            else if (ASSOCIATED(this%reactive_zone)) then                  !< If only reactive zone is associated
                this%concentrations(1)=1d-32                               !< Set free site as first concentration
                                                                            !! Default position when no mineral zone
            else
                error stop "Reactive zone not associated to solid chemistry" !< Error if no zone associations
            end if
        end subroutine
       
        
        !> Set volumetric fractions of minerals
        !! Physical Context: Volume fractions are dimensionless (0 ≤ φ ≤ 1) representing mineral volume / bulk volume
        !! Conservation: Σφᵢ + φ_pore = 1 where φ_pore is porosity
        !! Order: Kinetic minerals first, then equilibrium minerals
        subroutine set_vol_fracts(this,min_pore_vol_ratios)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: min_pore_vol_ratios(:)                 !< Volumetric fractions array (dimensionless)
            if (this%mineral_zone%num_minerals<size(min_pore_vol_ratios)) error stop "Dimension error in set_vol_fracts"
                                                                            !< Validate array size against number of minerals
            this%min_pore_vol_ratios=min_pore_vol_ratios
        end subroutine
        
        subroutine allocate_min_pore_vol_ratios(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%min_pore_vol_ratios)) then
                deallocate(this%min_pore_vol_ratios)
            end if
            allocate(this%min_pore_vol_ratios(this%mineral_zone%num_minerals))
            this%min_pore_vol_ratios=0d0
        end subroutine
        
        subroutine allocate_vol_fracts_microorgs(this)
            class(solid_chemistry_c) :: this
            allocate(this%vol_fracts_microorgs(this%biofilm%num_microorgs))
            this%vol_fracts_microorgs=0d0
        end subroutine
        
        subroutine set_react_surfaces(this,react_surfaces)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: react_surfaces(:)
            if (this%mineral_zone%num_minerals<size(react_surfaces)) error stop "Dimension error in set_react_surfaces"
            this%react_surfaces=react_surfaces
        end subroutine
        
        subroutine allocate_react_surfaces(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%react_surfaces)) then
                deallocate(this%react_surfaces)
            end if
            allocate(this%react_surfaces(this%mineral_zone%num_minerals))
            this%react_surfaces=0d0
        end subroutine
        
        subroutine allocate_conc_solids(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%concentrations)) then
                deallocate(this%concentrations)
            end if
            if (allocated(this%conc_old)) then
                deallocate(this%conc_old)
            end if
            if (allocated(this%conc_old_old)) then
                deallocate(this%conc_old_old)
            end if
            !if (associated(this%mineral_zone)) then
            !    allocate(this%concentrations(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else
                allocate(this%concentrations(this%num_solids),this%conc_old(this%num_solids),this%conc_old_old(this%num_solids))
            !end if
            !this%concentrations=0d0
        end subroutine
        
        !> Allocate memory for solid species activities array
        !! Physical Context: Activities of solid species (typically equal to concentrations for pure solids)
        !! For minerals: activity = 1 for pure phases, concentration for solid solutions
        !! For surface complexes: activity = concentration typically
        subroutine allocate_activities(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            if (allocated(this%activities)) then                           !< Check if already allocated
                deallocate(this%activities)                                !< Deallocate existing array
            end if
            !if (associated(this%mineral_zone)) then                       !< Commented mineral zone check
                !allocate(this%activities(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else                                                          !< End commented section
                allocate(this%activities(this%num_solids))                 !< Allocate activities array
            !end if                                                        !< End commented section
            this%activities=0d0                                            !< Initialize to zero
        end subroutine
        
        !> Allocate memory for logarithmic activity coefficients and their Jacobian matrix
        !! Mathematical Context: log_γᵢ and d(log_γᵢ)/d(log_cⱼ) for activity coefficient corrections
        !! For solids: usually log_γᵢ = 0 (ideal behavior) except for solid solutions
        !! Jacobian used in Newton-Raphson methods for non-ideal systems
        subroutine allocate_log_act_coeffs_solid_chem(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            if (allocated(this%log_act_coeffs)) then                       !< Check if already allocated
                deallocate(this%log_act_coeffs)                            !< Deallocate existing array
            end if
            if (allocated(this%log_Jacobian_act_coeffs)) then              !< Same for the Jacobian (allocated together below)
                deallocate(this%log_Jacobian_act_coeffs)
            end if
            !if (associated(this%mineral_zone)) then                       !< Commented mineral zone allocation
            !    allocate(this%log_act_coeffs(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else                                                          !< End commented section
                allocate(this%log_act_coeffs(this%num_solids),this%log_Jacobian_act_coeffs(this%num_solids,this%num_solids))
                                                                            !< Allocate both coefficient vector and Jacobian matrix
            !end if                                                        !< End commented section
            this%log_act_coeffs=0d0                                        !< Initialize coefficients to zero (ideal behavior)
            this%log_Jacobian_act_coeffs=0d0                               !< Initialize Jacobian to zero
        end subroutine
        
        !> Allocate memory for equivalents of exchangeable cations
        !! Physical Context: Equivalents [eq/L_bulk] for charge balance in cation exchange
        !! Dimension: Number of exchangeable cation categories (Ca²⁺, Mg²⁺, Na⁺, etc.)
        !! Conservation: Σ(equivalents) = CEC (Cation Exchange Capacity)
        subroutine allocate_equivalents(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            if (allocated(this%equivalents)) then                          !< Check if already allocated
                deallocate(this%equivalents)                               !< Deallocate existing array
            end if
            allocate(this%equivalents(this%reactive_zone%cat_exch_zone%num_exch_cats))
                                                                            !< Allocate based on number of exchangeable cation types
        end subroutine
        
        !> Associate solid chemistry object with a reactive zone
        !! Design Pattern: Pointer association for shared data access
        !! Validation: Ensures reactive zone has valid chemical system before association
        subroutine set_reactive_zone(this,reactive_zone)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            class(reactive_zone_c), intent(in), target :: reactive_zone    !< Reactive zone object to associate (target attribute for pointer)
            if (associated(reactive_zone%chem_syst)) then                  !< Validate reactive zone has chemical system
                this%reactive_zone=>reactive_zone                          !< Create pointer association
            else
                error stop "Reactive zone object is not associated to a chemical system"  !< Error if invalid zone
            end if
        end subroutine  
        
        !> Associate solid chemistry object with a mineral zone
        !! Design Pattern: Pointer association for mineral-specific data
        !! Validation: Ensures mineral zone has valid chemical system before association
        subroutine set_mineral_zone(this,mineral_zone)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            type(mineral_zone_c), intent(in), target :: mineral_zone       !< Mineral zone object to associate (target attribute for pointer)
            if (associated(mineral_zone%chem_syst)) then                   !< Validate mineral zone has chemical system
                this%mineral_zone=>mineral_zone
            else
                error stop "mineral zone object is not associated to a chemical system"
            end if
        end subroutine
        
       
        
        subroutine compute_activities_solids(this)
        !> Computes activities of solid species from concentrations and activity coefficients
            class(solid_chemistry_c) :: this
            this%activities=this%concentrations*(10**this%log_act_coeffs)
        end subroutine
        
       
        
        !> @brief Compute equivalents of adsorbed cations from concentrations
        !> @details Calculates the charge equivalents [eq/L_bulk] of adsorbed cations on exchange sites
        !> @par Mathematical Context
        !> equivalents_i = z_i × c_i where z_i is the valence and c_i is the concentration
        !> @par Physical Context
        !> - Equivalents represent the total exchangeable charge contributed by each cation type
        !> - Units: [eq/L_bulk] = [mol/L_bulk] × [charge/mol]
        !> - Conservation: Σ(equivalents_i) = CEC (Cation Exchange Capacity)
        !> @par Implementation
        !> - Loops through all exchangeable cation categories
        !> - Multiplies valence by concentration for each cation type
        !> - Assumes concentrations are stored after minerals and free site (position: num_minerals+1+i)
        !> @param[in,out] this Solid chemistry object (polymorphic class)
        !> @note Array indexing: concentrations(num_minerals+1+i) where i is cation index
        !> @see set_CEC, allocate_equivalents, compute_conc_surf_compl
        subroutine compute_equivalents(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (modified in-place)
            integer(kind=4) :: i                                           !< Loop index for exchangeable cation categories
            !> @par Algorithm
            !> For each exchangeable cation type i:
            !> 1. Get valence z_i from surface complex definition
            !> 2. Get concentration c_i from solid concentrations array
            !> 3. Calculate equivalents_i = z_i × c_i
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats           !< Loop over all exchangeable cation types
                this%equivalents(i)=this%reactive_zone%cat_exch_zone%surf_compl(1+i)%valence*this%concentrations(&
                this%reactive_zone%num_minerals+1+i)                       !< equivalents = valence × concentration
                                                                            !! Position: after minerals and free sites in concentrations array
            end do
        end subroutine
        
        !> Set Cation Exchange Capacity (CEC) with validation
        !! Physical Context: CEC [eq/L_bulk] represents total charge that can be balanced by exchangeable cations
        !! Typical values: 0.01-1.0 eq/L for soils and sediments
        !! Constraint: CEC ≥ 0 (cannot have negative charge capacity)
        subroutine set_CEC(this,CEC)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: CEC                                !< Cation Exchange Capacity [eq/L_bulk]
            if (CEC<0d0) then                                               !< Validate non-negative CEC
                error stop "CEC cannot be negative"                        !< Error for invalid CEC
            else
                this%CEC=CEC                                               !< Assign valid CEC value
            end if
        end subroutine
        
        !> Set activities of surface complexes including free sites
        !! Mathematical Context: Surface complexation equilibrium with site balance constraint
        !! Conservation: X_free + Σ(X_i-cation) = 1 where X represents site fractions
        !! Physical Context: Adsorbed cations compete for surface sites
        subroutine set_act_surf_compl(this,act_ads_cats)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(in) :: act_ads_cats(:)                    !< Activities of adsorbed cations (surface complexes)
            integer(kind=4) :: i                                           !< Loop index for exchangeable cations
            !> Set free site activity using site balance: X_free = 1 - Σ(X_cations)
            this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats)=1d0-SUM(act_ads_cats)
                                                                            !< Free site activity = 1 minus sum of occupied sites
            !> Set individual adsorbed cation activities
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats           !< Loop over exchangeable cation types
                this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats+i)=act_ads_cats(i)
                                                                            !< Assign activity to each adsorbed cation species
            end do
        end subroutine
        
        !> Set indices for solid species classification (variable vs constant activity)
        !! Programming Context: Separates solids into variable and constant activity categories
        !! Variable activity: Minerals that can dissolve/precipitate (equilibrium minerals)
        !! Constant activity: Pure phase minerals with fixed activity = 1
        subroutine set_indices_solids(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            integer(kind=4) :: i,j,k                                       !< Loop and counting indices
            j=0                                                             !< Counter for variable activity species
            k=0                                                             !< Counter for constant activity species
            !> Classify each mineral by its activity behavior
            do i=1,this%mineral_zone%num_minerals                          !< Loop over all minerals
                if (this%mineral_zone%chem_syst%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                                                                            !< Check if mineral has variable activity
                    j=j+1                                                   !< Increment variable activity counter
                    this%var_act_species_indices(j)=i                      !< Store index in variable activity array
                else                                                        !< If mineral has constant activity
                    k=k+1                                                   !< Increment constant activity counter
                    this%cst_act_species_indices(k)=i                      !< Store index in constant activity array
                end if
            end do
        end subroutine
        
    !> Update concentrations of adsorbed cations using Newton method corrections
    !! Mathematical Context: Newton-Raphson update: c_new = c_old + Δc where Δc is computed from Jacobian
    !! Physical Context: Surface complexation equilibrium solving with charge balance constraints
    !! Convergence: Iterative refinement until |Δc| < tolerance
        subroutine update_conc_ads_cats(this,conc_ads_cats,Delta_conc_ads_cats)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            real(kind=8), intent(inout) :: conc_ads_cats(:)                !< Concentration of adsorbed cations [mol/L_bulk]
            real(kind=8), intent(inout) :: Delta_conc_ads_cats(:)          !< Newton correction to adsorbed cation concentrations
    
            integer(kind=4) :: i                                           !< Loop index for cation types
            real(kind=8), allocatable :: conc_old(:)                       !< Backup of old concentrations for stability
    
            if (this%reactive_zone%CV_params%control_factor>1d0 .or. this%reactive_zone%CV_params%control_factor<0d0) then
                error stop "Control factor must be in (0,1)"
            end if
            conc_old=conc_ads_cats
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                if (conc_ads_cats(i)+Delta_conc_ads_cats(i)<=this%reactive_zone%CV_params%control_factor*conc_ads_cats(i)) then
                    conc_ads_cats(i)=this%reactive_zone%CV_params%control_factor*conc_ads_cats(i)
                else if (conc_ads_cats(i)+Delta_conc_ads_cats(i)>=conc_ads_cats(i)/this%reactive_zone%CV_params%control_factor) then
                    conc_ads_cats(i)=conc_ads_cats(i)/this%reactive_zone%CV_params%control_factor
                else
                    conc_ads_cats(i)=conc_ads_cats(i)+Delta_conc_ads_cats(i)
                end if
                Delta_conc_ads_cats(i)=conc_ads_cats(i)-conc_old(i)
            end do
        end subroutine
        
    !> Updates activity adsorbed cations in Newton method
        subroutine update_act_ads_cats(this,act_ads_cats,Delta_act_ads_cats)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(inout) :: act_ads_cats(:) !> concentration adsorbed cations
            real(kind=8), intent(inout) :: Delta_act_ads_cats(:) !> adsorbed cation concentration difference
    
            integer(kind=4) :: i
            real(kind=8), allocatable :: act_old(:)
    
            !> Validate control factor for Newton step size control
            if (this%reactive_zone%CV_params%control_factor>1d0 .or. this%reactive_zone%CV_params%control_factor<0d0) then
                error stop "Control factor must be in (0,1)"             !< Error if control factor outside valid range [0,1]
            end if
            act_old=act_ads_cats                                           !< Store old activities for backtracking if needed
            !> Newton step size control loop with backtracking for numerical stability
            do                                                              !< Infinite loop with internal exit condition
                !> Apply Newton correction with adaptive step size control
                do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats       !< Loop over all exchangeable cation types
                    !> Control Newton step size to prevent large oscillations or overshooting
                    if (act_ads_cats(i)+Delta_act_ads_cats(i)<=this%reactive_zone%CV_params%control_factor*act_ads_cats(i)) then
                        act_ads_cats(i)=this%reactive_zone%CV_params%control_factor*act_ads_cats(i)  !< Apply lower bound control
                    else if (act_ads_cats(i)+Delta_act_ads_cats(i)>=act_ads_cats(i)/this%reactive_zone%CV_params%control_factor)then
                        act_ads_cats(i)=act_ads_cats(i)/this%reactive_zone%CV_params%control_factor   !< Apply upper bound control
                    else
                        act_ads_cats(i)=act_ads_cats(i)+Delta_act_ads_cats(i)  !< Apply full Newton step if within bounds
                    end if
                    Delta_act_ads_cats(i)=act_ads_cats(i)-act_old(i)       !< Update correction based on actually applied step
                end do
                !> Check site balance constraint: Σ(X_i) < 1 (total site occupancy cannot exceed unity)
                if (SUM(act_ads_cats)>=1d0) then                           !< If total site occupancy violates physical constraint
                    act_ads_cats=act_old                                   !< Restore previous activities (backtrack)
                    Delta_act_ads_cats=Delta_act_ads_cats/2d0               !< Halve the correction step size
                    if (inf_norm_vec_real(Delta_act_ads_cats)<this%reactive_zone%CV_params%abs_tol) then
                        error stop "Delta_act_ads_cats is too small"       !< Error if correction becomes negligible
                    end if
                else                                                        !< If site balance constraint is satisfied
                    exit                                                    !< Accept Newton step and exit control loop
                end if
            end do
        end subroutine
        
        !> Compute surface complex concentrations from activities and activity coefficients
        !! Mathematical Context: c_i = a_i / γ_i where c=concentration, a=activity, γ=activity coefficient
        !! Physical Context: Converts thermodynamic activities to measurable concentrations
        !! For surface complexes: typically γ_i ≈ 1 (ideal behavior) unless strong electrostatic effects present
        subroutine compute_conc_surf_compl(this)
            !> Computes concentration of surface complexes from activities and activity coefficients
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            integer(kind=4) :: i                                           !< Loop index for surface complex species
            !> Convert activities to concentrations using activity coefficients
            do i=1,this%reactive_zone%cat_exch_zone%num_surf_compl          !< Loop over all surface complexes
                this%concentrations(this%num_solids-this%reactive_zone%cat_exch_zone%num_surf_compl+i)=&
                this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_surf_compl+i)/(&
                10**this%log_act_coeffs(this%num_solids-this%reactive_zone%cat_exch_zone%num_surf_compl+i))
                                                                            !< concentration = activity / 10^(log_activity_coefficient)
                                                                            !! Position: surface complexes stored at end of solid species array
            end do
            !this%activities(1)=this%concentrations(1)                     !< Commented free site assignment (temporary hack - "chapuza")
            !print *, this%concentrations                                  !< Commented debug output for concentrations
            !print *, this%activities                                      !< Commented debug output for activities
        end subroutine

        !> Compute total number of solid species in the system
        !! Count includes: reactive zone solids (equilibrium) + mineral zone kinetic minerals
        !! Used for array allocation and indexing throughout solid chemistry calculations
        subroutine compute_num_solids_solid_chem(this)
            class(solid_chemistry_c) :: this                               !< Solid chemistry object (polymorphic)
            this%num_solids=this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin
                                                                            !< Total = equilibrium solids + kinetic minerals
        end subroutine



!> Updates concentration solids in Newton method
subroutine update_conc_solids(this,Delta_c_s,control_factor)
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_s(:) !> solid concentration difference
    real(kind=8), intent(in) :: control_factor !> must \f$\in (0,1)\f$
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (control_factor>1d0 .or. control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,size(this%concentrations)
        if (this%concentrations(i)+Delta_c_s(i)<=control_factor*this%concentrations(i)) then
            this%concentrations(i)=control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_s(i)>=this%concentrations(i)/control_factor) then
            this%concentrations(i)=this%concentrations(i)/control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_s(i)
        end if
        Delta_c_s(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine

!> This subroutine checks for zero concentrations in a solid chemistry object
subroutine check_solid_chemistry(this,flag,indices)
    class(solid_chemistry_c) :: this
    !real(kind=8), intent(in) :: tolerance !> tolerance for concentrations of solids
    integer(kind=4), intent(out) :: flag !> 1 if no zero concentrations, 0 otherwise
    integer(kind=4), intent(out), allocatable :: indices(:) !> indices of zero concentrations
    
    integer(kind=4) :: num_new_ind_non_flow_species,i,j,k
    integer(kind=4), allocatable :: old_nf_ind(:),old_solid_ind(:)
    !type(species_c), allocatable :: new_ind_non_flow_species(:)
    type(solid_chemistry_c), allocatable :: new_solid_chems(:)
    real(kind=8), parameter :: epsilon=1d-9 !> arbitrario
    
    flag=1
    do i=1,this%reactive_zone%num_solids
        if (this%concentrations(i)<this%reactive_zone%CV_params%abs_tol) then
            flag=0
            indices=[indices,j] !> indices of zero concentration solids
        else
            continue
        end if
    end do
end subroutine

!> Computes concentration of minerals after a given time step
!! We assume minerals have constant activity
subroutine compute_conc_minerals_iter(this,Delta_t,lambda_r)
!> Arguments
    class(solid_chemistry_c) :: this !> solid chemistry object
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio; scales the eq amount to the aqueous contribution lambda_r*Re
!> Variables
    integer(kind=4) :: i !> counter minerals
    real(kind=8), parameter :: eps=1d-16 !> chapuza
!> Process
    do i=1,this%reactive_zone%num_minerals
        if (abs(this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i))<eps) then
            continue
        else
            this%concentrations(this%mineral_zone%num_minerals_kin+i)=this%concentrations(this%mineral_zone%num_minerals_kin+i)+&
                lambda_r*Delta_t*this%rxn_rate_eq(i)!/this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i)
            if (this%concentrations(this%mineral_zone%num_minerals_kin+i)<0d0) then
                this%concentrations(this%mineral_zone%num_minerals_kin+i)=0d0
            end if
        end if
    end do
    ! print *, "DEBUG compute_conc_minerals_iter: num_minerals_kin = ", this%mineral_zone%num_minerals_kin
    ! print *, "DEBUG: num_lin_kin_reacts   = ", this%mineral_zone%chem_syst%num_lin_kin_reacts
    ! print *, "DEBUG: num_redox_kin_reacts = ", this%mineral_zone%chem_syst%num_redox_kin_reacts
    ! print *, "DEBUG: num_kin_reacts       = ", this%mineral_zone%chem_syst%num_kin_reacts
    ! print *, "DEBUG: size(Sk)             = ", size(this%mineral_zone%chem_syst%Sk,1), " x ", &
    !     size(this%mineral_zone%chem_syst%Sk,2)
    ! print *, "DEBUG: size(rk_mean)        = ", size(this%rk_mean)
    ! print *, "DEBUG: rk_mean              = ", this%rk_mean
    ! print *, "DEBUG: Delta_t              = ", Delta_t
    do i=1,this%mineral_zone%num_minerals_kin
        ! print *, "DEBUG: kin mineral i=", i, " vol_fract=", this%min_pore_vol_ratios(i), &
        !     " ind_min_Sk=", this%mineral_zone%ind_min_Sk(i), " conc_before=", this%concentrations(i)
        if (abs(this%min_pore_vol_ratios(i))<eps) then
            !print *, "DEBUG: skipping i=", i, " (vol_fract < eps)"
            continue                                                            !< Skip to next iteration (no concentration update)
        else                                                                !< If mineral dissolution/precipitation is occurring
            !> Update mineral concentration using kinetic rate law
            ! print *, "DEBUG: Sk row range = ", &
            !     this%mineral_zone%chem_syst%num_lin_kin_reacts+this%mineral_zone%chem_syst%num_redox_kin_reacts+1, &
            !     " : ", this%mineral_zone%chem_syst%num_kin_reacts
            ! print *, "DEBUG: Sk column (ind_min_Sk) = ", this%mineral_zone%ind_min_Sk(i)
            ! print *, "DEBUG: Sk slice = ", this%mineral_zone%chem_syst%Sk(&
            !     this%mineral_zone%chem_syst%num_lin_kin_reacts+this%mineral_zone%chem_syst%num_redox_kin_reacts+1:&
            !     this%mineral_zone%chem_syst%num_kin_reacts, this%mineral_zone%ind_min_Sk(i))
            ! print *, "DEBUG: rk_mean slice = ", this%rk_mean(1:this%mineral_zone%num_minerals_kin)
            this%concentrations(i)=this%concentrations(i)+Delta_t*dot_product(this%mineral_zone%chem_syst%Sk(&
            this%mineral_zone%chem_syst%num_lin_kin_reacts+this%mineral_zone%chem_syst%num_redox_kin_reacts+1:&
            this%mineral_zone%chem_syst%num_kin_reacts,&
                this%mineral_zone%ind_min_Sk(i)),this%rxn_rate_kin_avg(1:this%mineral_zone%num_minerals_kin))!/&
            !this%min_pore_vol_ratios(i)
                                                                            !< dc/dt = (S^T × r) where S=stoichiometry, r=rates
            ! print *, "DEBUG: conc_after=", this%concentrations(i)
            if (this%concentrations(i)<0d0) then
                this%concentrations(i)=0d0                                   !< Prevent negative concentrations
                ! print *, "DEBUG: clamped to 0"
            end if
        end if
    end do
end subroutine

!> Modify mixing ratios for kinetic reaction rates to prevent negative concentrations
!! Mathematical Context: Adaptive time stepping for reactive transport to maintain positivity
!! Algorithm: If c + Δt×r < 0, reduce mixing ratio λ to ensure c + Δt×λ×r ≥ 0
!! Physical Context: Prevents unphysical negative concentrations during reactive mixing
subroutine modify_mix_ratios_reacts_solid_chem(this,mix_ratio_init,c_mix,Delta_t,r_tilde,mix_ratio_new,num_lump)
    !> This subroutine modifies the mixing ratios of the reaction amounts for a target
    !> AQUI DEBERIAS GUARDAR LOS NUEVOS LAMBDAS                           !< TODO: Should save new lambdas
        class(solid_chemistry_c), intent(in) :: this                       !< Solid chemistry object (input only)
        real(kind=8), intent(in) :: mix_ratio_init                         !< Initial mixing ratio of reaction rates
        real(kind=8), intent(in) :: c_mix(:)                             !< Concentration vector after transport
        real(kind=8), intent(in) :: Delta_t                                !< Time step size
        real(kind=8), intent(inout) :: r_tilde(:)                         !< reaction rate contributions
        real(kind=8), intent(out) :: mix_ratio_new                      !< New mixing ratio for kinetic reactions (output)
        integer(kind=4), intent(out) :: num_lump                           !< Number of lumping operations performed

        integer(kind=4) :: i                                               !< Loop index for species
        integer(kind=4) :: n_v                                            !< Number of variable activity species
        real(kind=8), parameter :: alpha=1.05                              !< Mixing ratio modification factor (5% increase)
        !real(kind=8), allocatable :: R_tilde(:)                            !< Reaction amount contributions after mixing [mol/L]

        !R_tilde=Delta_t*r_tilde
        n_v=this%reactive_zone%speciation_alg%num_var_act_species          !< Get number of variable activity species
        mix_ratio_new=mix_ratio_init                                 !< Initialize new mixing ratio (temporary - "chapuza")
        !flag=.false.                                                      !< Commented lumping flag initialization
        num_lump=0                                                          !< Initialize number of lumping operations
        i=1                                                                 !< Initialize species index
        !> Check for negative concentrations after mixing and adjust if necessary
        do                                                                  !< Infinite loop with exit condition
            !if (c_mix(i)+R_tilde_up(i)<=0d0) then                       !< Commented upstream check
                !error stop "Negative concentration after mixing upstream waters in subroutine reactive_mixing_iter_EI_kin_aq_anal_ideal_opt2"
            !if (c_mix(i)+R_tilde_up(i)+R_tilde_down(i)<=0d0 .and. R_tilde_down(i)<0d0) then  !< Commented combined check
            if (c_mix(i)+Delta_t*r_tilde(i)<=0d0) then                           !< Check if total concentration would become negative
                !flag=.true.                                               !< Commented lumping flag setting
                num_lump=num_lump+1                                        !< Increment lumping counter
                !print *, "Warning: negative concentration after mixing"   !< Commented warning output
                !print *, "Species index: ", i                             !< Commented species index output
                !R_tilde_down(i)=R_tilde_down(i)/2d0                       !< Commented simple halving (authentic hack)
                !R_tilde_down=R_tilde_down*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new)  !< Commented downstream reduction
                !R_tilde_up=R_tilde_up*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new)    !< Commented upstream reduction
                r_tilde=r_tilde*(1d0-alpha*mix_ratio_new)/(1d0-mix_ratio_new)     !< Reduce reaction contributions (hack)
                mix_ratio_new=mix_ratio_new*alpha                    !< Increase mixing ratio by factor α
                if (mix_ratio_new>1d0) then                     !< Check for excessive mixing ratio
                    !> Limit mixing ratio to prevent instability (hack)
                    mix_ratio_new=1d0                           !< Set maximum stable mixing ratio
                    !R_tilde_down=0d0                                      !< Commented downstream nullification
                    !R_tilde_up=0d0                                        !< Commented upstream nullification
                    r_tilde=0d0
                    exit
                end if
            else if (i<n_v) then
                i=i+1
            else
                !r_tilde=R_tilde/Delta_t !> we update kinetic reaction rate contributions
                exit
            end if
        end do
end subroutine

subroutine copy_solid_chemistry(this,solid_chemistry)
class(solid_chemistry_c) :: this !> solid chemistry object
type(solid_chemistry_c), intent(in) :: solid_chemistry !> solid chemistry object to assign
call this%copy_local_chemistry(solid_chemistry)
if (associated(solid_chemistry%reactive_zone)) then
    this%reactive_zone=>solid_chemistry%reactive_zone
else
    error stop "Reactive zone not associated with solid chemistry"
end if
if (associated(solid_chemistry%mineral_zone)) then
    this%mineral_zone=>solid_chemistry%mineral_zone
else
    error stop "Mineral zone not associated with solid chemistry"
end if
if (associated(solid_chemistry%tar)) then
    this%tar=>solid_chemistry%tar
else
    !> Source has no target yet (e.g. default solid chemistry copied before
    !> targets are linked via set_target). Nullify here; the spatial target
    !> is associated later in chemistry_m::set_target_solids_mesh.
    this%tar=>null()
end if
this%num_solids=solid_chemistry%num_solids
this%time=solid_chemistry%time
this%is_interface=solid_chemistry%is_interface
if (allocated(solid_chemistry%vol_fracts_microorgs)) then
    this%vol_fracts_microorgs=solid_chemistry%vol_fracts_microorgs
end if
if (allocated(solid_chemistry%min_pore_vol_ratios)) then
    this%min_pore_vol_ratios=solid_chemistry%min_pore_vol_ratios
end if
if (allocated(solid_chemistry%react_surfaces)) then
    this%react_surfaces=solid_chemistry%react_surfaces
end if
if (allocated(solid_chemistry%equivalents)) then
    this%equivalents=solid_chemistry%equivalents
end if
this%CEC=solid_chemistry%CEC
end subroutine

subroutine set_target(this,target_obj)
    class(solid_chemistry_c) :: this
    type(target_c), intent(in), target :: target_obj
    this%tar=>target_obj
end subroutine

subroutine allocate_rxn_rates_and_amounts_solid_chem(this)
    class(solid_chemistry_c) :: this !> solid chemistry object
    integer(kind=4) :: n_eq !> number of solid equilibrium reactions (minerals + exchange)
    if (.not. allocated(this%rxn_rate_eq)) then
        if (.not. associated(this%reactive_zone)) then
            error stop "Reactive zone not associated with solid chemistry"
        else
            n_eq=this%reactive_zone%num_minerals+this%reactive_zone%cat_exch_zone%num_exch_cats
            !> full eq-rate history: aqueous_chemistry_c%allocate_rxn_rates_and_amounts skips this block
            !> once rxn_rate_eq is allocated, so the *_old arrays must be allocated here too
            allocate(this%rxn_rate_eq(n_eq))
            allocate(this%rxn_amnt_eq(n_eq))
            allocate(this%cum_rxn_amnt_eq(n_eq))
            allocate(this%rxn_rate_eq_old(n_eq))
            allocate(this%rxn_rate_eq_old_old(n_eq))
            allocate(this%rxn_rate_eq_old_old_old(n_eq))
            this%rxn_rate_eq=0d0 !> by default
            this%rxn_amnt_eq=0d0 !> by default
            this%cum_rxn_amnt_eq=0d0 !> by default
            this%rxn_rate_eq_old=0d0 !> by default
            this%rxn_rate_eq_old_old=0d0 !> by default
            this%rxn_rate_eq_old_old_old=0d0 !> by default
        end if
    end if
    if (.not. allocated(this%rxn_amnt_kin)) then
        if (.not. associated(this%mineral_zone)) then
            error stop "Mineral zone not associated with solid chemistry"
        else
            allocate(this%rxn_amnt_kin(this%mineral_zone%num_minerals_kin))
            allocate(this%cum_rxn_amnt_kin(this%mineral_zone%num_minerals_kin))
            allocate(this%rxn_rate_kin(this%mineral_zone%num_minerals_kin))
            allocate(this%rxn_rate_kin_old(this%mineral_zone%num_minerals_kin))
            allocate(this%rk_old_old(this%mineral_zone%num_minerals_kin))
            allocate(this%rk_old_old_old(this%mineral_zone%num_minerals_kin))
            allocate(this%rxn_rate_kin_avg(this%mineral_zone%num_minerals_kin))
            this%rxn_amnt_kin=0d0 !> by default
            this%cum_rxn_amnt_kin=0d0 !> by default
            this%rxn_rate_kin=0d0 !> by default
            this%rxn_rate_kin_old=0d0 !> by default
            this%rk_old_old=0d0 !> by default
            this%rk_old_old_old=0d0 !> by default
            this%rxn_rate_kin_avg=0d0 !> by default
        end if
    end if
end subroutine

!> \brief Deallocates the equilibrium and kinetic reaction rate and amount arrays (the ones
!>        allocate_rxn_rates_and_amounts_solid_chem creates), so that a following call to it lays
!>        them out again for the CURRENT reactive and mineral zones (it skips arrays already allocated).
subroutine deallocate_rxn_rates_and_amounts_solid_chem(this)
    class(solid_chemistry_c) :: this !> solid chemistry object
    if (allocated(this%rxn_rate_eq)) deallocate(this%rxn_rate_eq)
    if (allocated(this%rxn_amnt_eq)) deallocate(this%rxn_amnt_eq)
    if (allocated(this%cum_rxn_amnt_eq)) deallocate(this%cum_rxn_amnt_eq)
    if (allocated(this%rxn_rate_eq_old)) deallocate(this%rxn_rate_eq_old)
    if (allocated(this%rxn_rate_eq_old_old)) deallocate(this%rxn_rate_eq_old_old)
    if (allocated(this%rxn_rate_eq_old_old_old)) deallocate(this%rxn_rate_eq_old_old_old)
    if (allocated(this%rxn_amnt_kin)) deallocate(this%rxn_amnt_kin)
    if (allocated(this%cum_rxn_amnt_kin)) deallocate(this%cum_rxn_amnt_kin)
    if (allocated(this%rxn_rate_kin)) deallocate(this%rxn_rate_kin)
    if (allocated(this%rxn_rate_kin_old)) deallocate(this%rxn_rate_kin_old)
    if (allocated(this%rk_old_old)) deallocate(this%rk_old_old)
    if (allocated(this%rk_old_old_old)) deallocate(this%rk_old_old_old)
    if (allocated(this%rxn_rate_kin_avg)) deallocate(this%rxn_rate_kin_avg)
end subroutine

!> \brief (Re)initialises a solid chemistry object from its associated reactive zone and mineral
!>        zone: computes the number of solids and allocates every state array from scratch.
!> \details Allocates the solid-species arrays (concentrations and their two previous-step copies,
!>          activities, log activity coefficients and their Jacobian), the exchanger equivalents, the
!>          mineral volume ratios and reactive surfaces, and the equilibrium and kinetic reaction
!>          rate and amount arrays. Solid layout: the minerals (equilibrium ones of the reactive zone,
!>          then kinetic ones of the mineral zone) lead, the surface complexes trail. Minerals are
!>          pure phases (c = a = 1, log gamma = 0); everything else starts at zero, CEC included (the
!>          exchanger state cannot be derived from the zones alone: set it afterwards if needed).
!>          Idempotent: already allocated arrays are released first, so the object can be re-pointed
!>          to other zones and re-initialised any number of times (e.g. the persistent scratch solid
!>          chemistry of a zone water, re-initialised at every within-zone mix); previous state is lost.
!> \pre reactive_zone and mineral_zone associated (set_reactive_zone, set_mineral_zone)
subroutine initialise_solid_chemistry(this)
    class(solid_chemistry_c) :: this !> solid chemistry object
    integer(kind=4) :: n_min !> number of minerals: equilibrium (reactive zone) + kinetic (mineral zone)
    if (.not. associated(this%reactive_zone)) error stop "initialise_solid_chemistry: reactive zone not associated"
    if (.not. associated(this%mineral_zone)) error stop "initialise_solid_chemistry: mineral zone not associated"
    call this%compute_num_solids_solid_chem()
    !> the array allocators below release an already allocated array themselves; the reaction
    !> rate/amount allocator keeps existing arrays, so they are released here to follow the zones
    call this%deallocate_rxn_rates_and_amounts_solid_chem()
    call this%allocate_conc_solids()
    call this%allocate_activities()
    call this%allocate_log_act_coeffs_solid_chem()
    call this%allocate_equivalents()
    call this%allocate_min_pore_vol_ratios()
    call this%allocate_react_surfaces()
    call this%allocate_rxn_rates_and_amounts_solid_chem()
    n_min=this%reactive_zone%num_minerals+this%mineral_zone%num_minerals_kin
    this%concentrations=0d0
    this%concentrations(1:n_min)=1d0
    this%conc_old=this%concentrations
    this%conc_old_old=this%concentrations
    this%activities(1:n_min)=1d0
    this%equivalents=0d0
    this%CEC=0d0
end subroutine

function get_tar_id(this) result(tar_id)
    class(solid_chemistry_c) :: this
    integer(kind=4) :: tar_id
    if (associated(this%tar)) then
        tar_id=this%tar%id
    else
        error stop "Target not associated with solid chemistry"
    end if
end function

subroutine set_time(this, time_val)
    class(solid_chemistry_c) :: this
    real(kind=8), intent(in) :: time_val
    if (time_val < 0d0) then
        error stop "Time cannot be negative"
    end if
    this%time = time_val
end subroutine

!=========================================================================
! BIOFILM: KOZENY-CARMAN PERMEABILITY REDUCTION (bioclogging)
!=========================================================================
!> \brief Compute relative permeability reduction due to bioclogging
!> \details Generalized Kozeny-Carman: k_b/k_0 = (n_b / n_0)^m, where
!> n_b = max(porosity_biofilm, biofilm%porosity_min) (local state / parameter),
!> n_0 = porosity_initial (local state) and m = biofilm%KC_exponent (parameter).
!> Requires the biofilm parameter object to be associated.
!> \param[in] this Solid chemistry object (carries the biofilm local state)
!> \return k_ratio Permeability ratio k_b/k_0 [-]
!> Ref: Hommel et al. (2018) Transport in Porous Media; Thullner et al. (2002) WRR
pure real(kind=8) function compute_permeability_reduction(this) result(k_ratio)
    class(solid_chemistry_c), intent(in) :: this                            !< Solid chemistry object
    real(kind=8) :: n_b                                                     !< Effective biofilm porosity

    n_b = max(this%porosity_biofilm, this%biofilm%porosity_min)

    if (this%porosity_initial > 0d0) then
        k_ratio = (n_b / this%porosity_initial) ** this%biofilm%KC_exponent
    else
        k_ratio = 1d0
    end if
end function

!=========================================================================
! BIOFILM: DETACHMENT RATE
!=========================================================================
!> \brief Compute total biofilm detachment rate [m/s]
!> \details Combines continuous erosion and sloughing:
!>   r_det = k_det * Lf + k_sl * max(0, Lf - Lf_max) + k_det_shear * tau^gamma * Lf
!> Thickness Lf = this%thickness (local state); the coefficients k_det, thickness_max,
!> k_sloughing, k_det_shear and det_shear_exponent are biofilm parameters.
!> \param[in] this Solid chemistry object (carries the biofilm local state)
!> \param[in] tau Wall shear stress [Pa] (optional; for shear-dependent term)
!> \return r_det Total detachment rate [m/s]
!> Ref: Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314; Abbas et al. (2012) MBE 9:215
pure real(kind=8) function compute_detachment_rate(this, tau) result(r_det)
    class(solid_chemistry_c), intent(in) :: this                            !< Solid chemistry object
    real(kind=8), intent(in), optional :: tau                               !< Wall shear stress [Pa]

    ! Continuous erosion: r = k_det * Lf
    r_det = this%biofilm%k_det * this%thickness

    ! Sloughing: r = k_sl * max(0, Lf - Lf_max)
    if (this%thickness > this%biofilm%thickness_max .and. this%biofilm%k_sloughing > 0d0) then
        r_det = r_det + this%biofilm%k_sloughing * (this%thickness - this%biofilm%thickness_max)
    end if

    ! Shear-dependent: r = delta * tau^gamma * Lf
    if (present(tau) .and. this%biofilm%k_det_shear > 0d0) then
        r_det = r_det + this%biofilm%k_det_shear * (tau ** this%biofilm%det_shear_exponent) * this%thickness
    end if
end function

!=========================================================================
! BIOFILM: THIELE MODULUS
!=========================================================================
!> \brief Compute the biofilm Thiele modulus (phi^2)
!> \details phi^2 = mu_max * X * Lf^2 / (D_eff * S_bulk * Y), with Lf = this%thickness
!> (local state). Low phi^2 => reaction-limited (compact); high phi^2 => diffusion-limited.
!> \param[in] this Solid chemistry object (carries the biofilm local state)
!> \param[in] mu_max Maximum specific growth rate [1/s]
!> \param[in] X Biomass density [M/L^3]
!> \param[in] D_eff Effective diffusion coefficient [m^2/s]
!> \param[in] S_bulk Bulk substrate concentration [M/L^3]
!> \param[in] Y Yield coefficient [-]
!> \return phi2 Thiele modulus phi^2 [-]
!> Ref: Picioreanu et al. (1998) Biotechnol. Bioeng. 57:718
pure real(kind=8) function compute_thiele_modulus(this, mu_max, X, D_eff, S_bulk, Y) result(phi2)
    class(solid_chemistry_c), intent(in) :: this                            !< Solid chemistry object
    real(kind=8), intent(in) :: mu_max                                      !< Max specific growth rate [1/s]
    real(kind=8), intent(in) :: X                                           !< Biomass density [M/L^3]
    real(kind=8), intent(in) :: D_eff                                       !< Effective diffusion [m^2/s]
    real(kind=8), intent(in) :: S_bulk                                      !< Bulk substrate concentration [M/L^3]
    real(kind=8), intent(in) :: Y                                           !< Yield coefficient [-]

    if (D_eff > 0d0 .and. S_bulk > 0d0 .and. Y > 0d0) then
        phi2 = mu_max * X * this%thickness**2 / (D_eff * S_bulk * Y)
    else
        phi2 = 0d0
    end if
end function

!=========================================================================
! SET IS INTERFACE
!=========================================================================
!> \brief Set the interface flag for the solid chemistry object
!> \param[inout] this Solid chemistry object
subroutine set_is_interface(this,flag)
    class(solid_chemistry_c), intent(inout) :: this
    logical, intent(in) :: flag
    this%is_interface = flag
end subroutine set_is_interface



end module