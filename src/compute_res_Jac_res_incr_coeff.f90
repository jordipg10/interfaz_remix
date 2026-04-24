!> \file compute_res_Jac_res_incr_coeff.f90
!> \brief Computes residual and its Jacobian using incremental coefficients (finite difference method)
!> \details This subroutine calculates both the residual vector and its Jacobian matrix for the Newton-Raphson
!> iteration used in solving chemical equilibrium problems. The Jacobian is computed numerically using finite
!> differences (incremental coefficients method) rather than analytically.
!>
!> The incremental coefficient method perturbs each primary species concentration by a small amount (eps),
!> recomputes the residual, and uses the finite difference to approximate the Jacobian:
!> \f[
!> \frac{\partial r_i}{\partial c_{1,j}} \approx \frac{r_i(c_1 + \epsilon \cdot e_j) - r_i(c_1)}{\epsilon}
!> \f]
!> where:
!> - \f$r_i\f$ is the i-th residual equation (mass balance or constraint) [-]
!> - \f$c_{1,j}\f$ is the j-th primary species concentration [M/L³]
!> - \f$\epsilon\f$ is a small perturbation parameter (eps) [-]
!> - \f$e_j\f$ is the j-th unit vector [-]
!>
!> After computing the Jacobian, the subroutine restores the unperturbed state and recomputes secondary
!> quantities (ionic strength, activity coefficients, activities) for consistency.
!>
!> \param[in,out] this Aqueous chemistry object containing speciation algebra, reactive zone, and chemical system data (state updated during computation) [-]
!> \param[in] c2 Secondary species concentrations for unperturbed residual calculation (CHAPUZA - workaround, dimension = num_eq_reactions) [M/L³]
!> \param[in] indices_icon Object containing indices of different icon (initial condition) types (e.g., minerals, gases, aqueous species) [-]
!> \param[in] n_icon Number of each icon type (array with counts for different categories) [-]
!> \param[in] indices_constrains Indices of constrained species or equations (for mass balance or other constraints) [-]
!> \param[in] ctot Total concentrations for constrained species (used in residual calculation) [M/L³]
!> \param[out] res Residual vector in Newton-Raphson iteration (difference between computed and target values) [-]
!> \param[out] Jac_res Jacobian matrix of residual with respect to primary species concentrations ∂r/∂c₁ [-]

subroutine compute_res_Jac_res_incr_coeff(this,c2,indices_icon,n_icon,indices_constrains,ctot,res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    use arrays_m, only: int_array_c !> Import aqueous chemistry class and integer array class for storing indices
    implicit none !> Enforce explicit variable declarations
    !> Pre-process
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone, speciation algebra, concentrations, and activities [-]
    real(kind=8), intent(in) :: c2(:) !> Secondary species concentrations for unperturbed residual (CHAPUZA - workaround, dimension = num_eq_reactions) [M/L³]
    class(int_array_c), intent(in) :: indices_icon !> Object containing indices of different icon (initial condition) types [-]
    integer(kind=4), intent(in) :: n_icon(:) !> Number of each icon type (array with counts for minerals, gases, aqueous species, etc.) [-]
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> Indices of constrained species or equations (2D array for complex constraint mapping) [-]
    real(kind=8), intent(in) :: ctot(:) !> Total concentrations for constrained species (target values for mass balance constraints) [M/L³]
    real(kind=8), intent(out) :: res(:) !> Residual vector in Newton-Raphson iteration (unperturbed residual) [-]
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian matrix of residual: ∂r/∂c₁ (computed using finite differences) [-]
    
    real(kind=8), allocatable :: c1_pert(:) !> Perturbed primary species concentrations (c₁ + ε·eⱼ) for finite difference calculation [M/L³]
    real(kind=8), allocatable :: res_pert(:) !> Perturbed residual vector (residual computed with perturbed primary concentrations) [-]
    real(kind=8), allocatable :: c1(:) !> Unperturbed primary species concentrations (stored for restoration after perturbations) [M/L³]
    real(kind=8), allocatable :: c2_pert(:) !> Perturbed secondary species concentrations (recomputed from perturbed c₁) [M/L³]
    real(kind=8), allocatable :: log_c2k(:) !> Logarithm of secondary species concentrations at iteration k (not used in current implementation) [-]
    real(kind=8), allocatable :: log_c2(:) !> Logarithm of secondary species concentrations (not used in current implementation) [-]
    real(kind=8), allocatable :: dc2_dc1(:,:) !> Jacobian of secondary concentrations with respect to primary concentrations ∂c₂/∂c₁ (not used in current implementation) [-]
    real(kind=8), allocatable :: out_prod(:,:) !> Outer product matrix for potential analytical calculations (not used in current implementation) [-]
    real(kind=8), allocatable :: Delta_c1(:) !> Newton step: Δc₁ = c₁^(i+1) - c₁^i (not used in current implementation) [M/L³]
    real(kind=8), allocatable :: abs_tol_res(:) !> Absolute tolerances for residuals in Newton-Raphson (not used in current implementation) [-]
    real(kind=8), allocatable :: mat_lin_syst(:,:) !> Matrix for linear system in potential analytical method (not used in current implementation) [-]
    real(kind=8), allocatable :: Se_aq_comp(:,:) !> Stoichiometric matrix for aqueous equilibrium reactions (not used in current implementation) [-]
    real(kind=8), allocatable :: K(:) !> Equilibrium constants for reactions (not used in current implementation) [-]
    real(kind=8), allocatable :: u_aq(:) !> Aqueous activities or log activities (not used in current implementation) [-]
    real(kind=8), allocatable :: z2(:) !> Charges of secondary species (not used in current implementation) [-]
    integer(kind=4) :: i !> Loop counter for residual equations (not used in current implementation) [-]
    integer(kind=4) :: j !> Loop counter for primary species (used for Jacobian column iteration) [-]
    integer(kind=4) :: ind_eqn !> Index of current equation in residual system (not used in current implementation) [-]
    integer(kind=4) :: niter_Picard !> Number of Picard iterations (not used in current implementation) [-]
    integer(kind=4) :: ind_cstr !> Index of constraint equation (not used in current implementation) [-]
    integer(kind=4) :: niter !> Number of iterations for Picard method when computing perturbed c₂ [-]
    integer(kind=4), allocatable :: indices_Jac(:) !> Vector of icon indices extracted from indices_icon object [-]
    integer(kind=4), allocatable :: cols(:) !> Column indices (not used in current implementation) [-]
    integer(kind=4), allocatable :: ind_aq_species(:) !> Indices of aqueous species (not used in current implementation) [-]
    integer(kind=4), allocatable :: counters(:) !> Counters for various species types (not used in current implementation) [-]
    logical :: flag_gas !> Flag indicating presence of gas phase (not used in current implementation) [-]
    logical :: flag_min !> Flag indicating presence of minerals (not used in current implementation) [-]
    logical :: flag_wat !> Flag indicating presence of water (not used in current implementation) [-]
    logical :: CV_flag !> Convergence flag from Picard iteration for perturbed secondary concentrations [-]
    
    allocate(c1_pert(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),res_pert(& !> Allocate perturbed primary concentration array (size = number of primary species)
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> Allocate perturbed residual array (size = number of primary species)
    
    indices_Jac=indices_icon%get_vector_int() !> Extract vector of icon indices from indices_icon object for use in residual computation [-]
    call this%compute_res_init(indices_icon,n_icon,indices_constrains,ctot,res) !> Compute initial (unperturbed) residual vector using current primary concentrations and constraints [-]
    c1=this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> Store unperturbed primary species concentrations for later restoration [M/L³]
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Loop over all primary species to compute each column of Jacobian
        !> Compute perturbation: perturb j-th primary species concentration
        c1_pert=c1 !> Initialize perturbed concentration vector to unperturbed values [M/L³]
        c1_pert(j)=c1_pert(j)+this%solid_chemistry%reactive_zone%CV_params%eps !> Add small perturbation ε to j-th primary species concentration [M/L³]
        !> Set perturbed primary concentrations in chemistry object
        call this%set_conc_prim_species(c1_pert) !> Update primary species concentrations in aqueous chemistry object to perturbed values [M/L³]
        !> Compute perturbed secondary concentrations from perturbed primary concentrations
        call this%compute_c2_from_c1_Picard(c1_pert,c2,c2_pert,niter,CV_flag) !> Use Picard iteration to compute secondary aqueous species concentrations from perturbed c₁ (initial guess from c2) [M/L³]
        !call this%compute_activities() !> Commented out: computation of all activities (aqueous, minerals, gases) from concentrations
        !call this%compute_pH() !> Commented out: computation of pH = -log₁₀(aH⁺) from hydrogen ion activity
        !call this%compute_salinity() !> Commented out: computation of salinity from ionic composition
        !call this%compute_alkalinity() !> Commented out: computation of alkalinity from carbonate species and other contributors
    !> Compute perturbed residual with perturbed concentrations
        call this%compute_res_init(indices_icon,n_icon,indices_constrains,ctot,res_pert) !> Compute residual vector with perturbed primary concentrations and updated secondary concentrations [-]
        Jac_res(:,j)=(res_pert-res)/this%solid_chemistry%reactive_zone%CV_params%eps !> Compute j-th column of Jacobian using finite difference: ∂r/∂c₁ⱼ ≈ (r_pert - r)/ε [-]
    end do !> End loop over primary species
    call this%set_conc_prim_species(c1) !> Restore unperturbed primary species concentrations in aqueous chemistry object [M/L³]
    call this%set_conc_sec_aq_species(c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)) !> Restore unperturbed secondary aqueous species concentrations (extract aqueous species from full c2 array) [M/L³]
    call this%compute_ionic_strength() !> Recompute ionic strength I = 0.5·Σ(cᵢ·zᵢ²) using restored unperturbed concentrations [M/L³]
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_strength,this%params_aq_sol,this%log_act_coeffs) !> Recompute log₁₀(activity coefficients) for aqueous species using restored ionic strength and activity coefficient model parameters [-]
    call this%compute_activities_aq() !> Recompute aqueous species activities aᵢ = γᵢ·cᵢ using restored concentrations and activity coefficients (CHAPUZA - workaround) [M/L³]
    call this%compute_log_act_coeff_wat() !> Recompute log₁₀(activity coefficient) for water using restored state [-]
end subroutine !> End of compute_res_Jac_res_incr_coeff subroutine