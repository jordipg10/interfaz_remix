!> \file compute_conc_surf_ideal_bin.f90
!> \brief Computes surface complex activities from cation concentrations using analytical solution for ideal two-cation exchange
!> \details This subroutine computes adsorbed cation activities (site fractions) for a two-cation ion exchange system
!> using an analytical (explicit) solution. Unlike the Newton-Raphson approach, this provides a direct calculation
!> based on equilibrium constants and aqueous cation concentrations.
!>
!> For a two-cation exchange system with cations Ca²⁺ and Mg²⁺ (or similar), the analytical solution is:
!> \f[
!> a_2^{ads} = \frac{c_2 \cdot K_1}{K_1 K_2 + c_1 K_2 + c_2 K_1}
!> \f]
!> \f[
!> a_1^{ads} = \frac{c_1 K_2 \cdot a_2^{ads}}{K_1 c_2}
!> \f]
!> \f[
!> a_0^{ads} = 1 - (a_1^{ads} + a_2^{ads})
!> \f]
!> where:
!> - \f$ a_i^{ads} \f$ = activity (site fraction) of adsorbed cation i [-]
!> - \f$ c_i \f$ = aqueous concentration of cation i [C]
!> - \f$ K_i \f$ = inverse equilibrium constant (1/K_eq) for cation exchange half-reaction i [-]
!> - \f$ a_0^{ads} \f$ = activity of free (unoccupied) exchange sites [-]
!>
!> Assumptions:
!> - Ideal activity coefficients for aqueous cations: γᵢ = 1
!> - Exactly two cation exchange half-reactions
!> - Mass balance constraint: a₀ + a₁ + a₂ = 1
!>
!> \param[in,out] this Solid chemistry object containing reactive zone, cation exchange zone, and equilibrium reactions
!> \param[in] conc_cats Aqueous cation concentrations (dimension = 2 for two cations) [C]

subroutine compute_conc_surf_ideal_bin(this,conc_cats)
    use solid_chemistry_m, only: solid_chemistry_c !> Import solid chemistry class for reactive zone and surface complexation
    implicit none !> Enforce explicit variable declarations
    class(solid_chemistry_c) :: this !> Solid chemistry object containing reactive zone, cation exchange zone, and chemical system [-]
    real(kind=8), intent(in) :: conc_cats(:) !> Aqueous cation concentrations (dimension = number of cation exchange half-reactions, typically 2) [C]
    !real(kind=8), intent(in) :: act_ads_cats_ig(:) !> Commented: Initial guess for surface complex activities (not needed for analytical solution)
    !integer(kind=4), intent(out) :: niter !> Commented: Number of iterations (not applicable for analytical solution)
    !logical, intent(out) :: CV_flag !> Commented: Convergence flag (not applicable for analytical solution)
    
    integer(kind=4) :: n_ads_cats !> Number of adsorbed cations (not used in current implementation) [-]
    integer(kind=4) :: n_eq !> Total number of equilibrium reactions in system [-]
    integer(kind=4) :: n_sec_aq !> Number of secondary aqueous species (not used in this routine) [-]
    integer(kind=4) :: n_gas_eq !> Number of gas equilibrium reactions [-]
    integer(kind=4) :: i !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: j !> Loop counter (not used in current implementation) [-]
    real(kind=8) :: K_1 !> Inverse equilibrium constant for first cation exchange half-reaction: K₁ = 1/K_eq,1 [-]
    real(kind=8) :: K_2 !> Inverse equilibrium constant for second cation exchange half-reaction: K₂ = 1/K_eq,2 [-]
    real(kind=8), allocatable :: act_ads_cats_old(:) !> Old adsorbed cation activities (not used in current implementation) [-]
    real(kind=8), allocatable :: log_act_ads_cats_new(:) !> Logarithm of new adsorbed cation activities (not used in current implementation) [-]
    real(kind=8), allocatable :: act_ads_cats_new(:) !> New adsorbed cation activities (not used in current implementation) [-]
!> Pre-processing
    !n_p_aq=this%reactive_zone%speciation_alg%num_aq_prim_species !> Commented: Number of aqueous primary species (not used)
    n_eq=this%reactive_zone%speciation_alg%num_eq_reactions !> Extract total number of equilibrium reactions from speciation algebra object
    n_gas_eq=this%reactive_zone%gas_phase%num_gases_eq !> Extract number of gas equilibrium reactions from gas phase object
    !n_ads_cats=this%reactive_zone%cat_exch_zone%num_exch_cats !> Commented: Number of exchange cations (not used)
    !allocate(Delta_act_ads_cats(n_ads_cats),residual(n_ads_cats),Jacobian(n_ads_cats,n_ads_cats)) !> Commented: Arrays for Newton iteration (not needed)
    !allocate(log_act_ads_cats_new(n_ads_cats)) !> Commented: Array for logarithmic activities (not needed)
    K_1 = 1d0/this%reactive_zone%chem_syst%eq_reacts(& !> Extract equilibrium constant for first cation exchange half-reaction
        this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq-1))%eq_cst !> Index: second-to-last non-gas equilibrium reaction, then compute inverse: K₁ = 1/K_eq,1
    K_2 = 1d0/this%reactive_zone%chem_syst%eq_reacts(& !> Extract equilibrium constant for second cation exchange half-reaction
        this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq))%eq_cst !> Index: last non-gas equilibrium reaction, then compute inverse: K₂ = 1/K_eq,2
!> Process: Analytical solution for two-cation ion exchange equilibrium
    this%activities(this%num_solids)=conc_cats(2)*K_1/(K_1*K_2+conc_cats(1)*K_2+conc_cats(2)*K_1) !> Compute activity of second adsorbed cation: a₂ = c₂·K₁/(K₁K₂ + c₁K₂ + c₂K₁)
    this%activities(this%num_solids-1)=conc_cats(1)*K_2*this%activities(this%num_solids)/(K_1*conc_cats(2)) !> Compute activity of first adsorbed cation: a₁ = (c₁·K₂·a₂)/(K₁·c₂)
    this%activities(this%num_solids-2)=1d0-SUM(this%activities(this%num_solids-1:this%num_solids)) !> Compute activity of free (unoccupied) exchange sites: a₀ = 1 - (a₁ + a₂) to satisfy mass balance
!> Post-processing
    !call this%set_act_surf_compl(act_ads_cats_new) !> Commented: Store surface complex activities (already stored directly in this%activities)
    call this%compute_conc_surf_compl() !> Compute surface complex concentrations from activities (convert site fractions to concentrations)
 end subroutine