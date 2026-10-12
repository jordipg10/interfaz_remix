!> \file compute_c2_from_c1_Picard.f90
!> \brief Computes secondary species concentrations from primary species using Picard iteration
!> \details This subroutine solves the chemical equilibrium speciation problem using the
!> Picard fixed-point iteration method combined with the mass action law. The method is
!> suitable for systems with both aqueous and solid primary species.
!>
!> **Mathematical Formulation**:
!>
!> The mass action law for secondary species in matrix form:
!> \f[
!> \log c_2 = S_{e,1}^* (\log \gamma_1 + \log c_1) + \log \tilde{K} - \log \gamma_2
!> \f]
!> where:
!> - \f$c_2\f$: Secondary species concentrations [M]
!> - \f$c_1\f$: Primary species concentrations [M]
!> - \f$\gamma_1, \gamma_2\f$: Activity coefficients (primary, secondary) [-]
!> - \f$S_{e,1}^*\f$: Modified stoichiometric matrix [n_sec × n_prim]
!> - \f$\tilde{K}\f$: Modified equilibrium constants [-]
!>
!> **Picard Iteration Algorithm** (optimised as compute_c2v_from_c1_Picard):
!>
!> Only the activity coefficients of the AQUEOUS species change along the iteration (they depend on
!> the ionic strength, i.e. on the aqueous secondary concentrations themselves); \f$c_1\f$,
!> \f$\log\tilde{K}\f$ and the activity coefficients of the non-aqueous primary (free site, kinetic
!> solids and gases) and secondary (minerals, surface complexes, gases) species are fixed, so their
!> whole contribution is assembled ONCE before the loop (log_c2_cst). The iteration works on local
!> arrays only (no object updates, allocations or array temporaries per iteration):
!> 1. Ionic strength \f$I = I_{fixed} + \tfrac12\sum_i z_i^2 c_{2,i}\f$ over the aqueous secondary
!>    species, with \f$I_{fixed}\f$ (aqueous primary species, from the object's concentrations) and the
!>    squared valences precomputed
!> 2. \f$\log_{10}\gamma\f$ of the dissolved solids from \f$I\f$ with the per-species model parameters
!>    gathered once: \f$\log\gamma_i = -\alpha_i A z_i^2\sqrt{I}/(1+\beta_i\sqrt{I}) + \gamma_i I\f$
!>    (the formula of aq_phase_c%compute_log_act_coeffs_aq_phase; the water's \f$\log\gamma\f$ is fixed)
!> 3. Aqueous \f$\gamma\f$ corrections of the mass action law, accumulated column by column over the
!>    aqueous primary columns of \f$S_{e,1}^*\f$, then the exponentiation
!> 4. Convergence check: \f$\|\frac{c_2^{(k+1)} - c_2^{(k)}}{c_2^{(k)}}\|_\infty < \text{rel\_tol}\f$ AND
!>    \f$\|\frac{\log c_2^{(k+1)} - \log c_2^{(k)}}{\log c_2^{(k)}}\|_\infty < 10^{\text{log\_rel\_tol}}\f$
!>
!> The object is updated once after convergence: ionic strength and aqueous log activity
!> coefficients at the converged \f$c_2\f$, the secondary concentrations (aqueous, solid, gas),
!> activities, the water activity coefficient, the salinity and the gas partial pressures.
!>
!> **Differences with the previous version** (both fixes make the result independent of the species
!> order of the water, as compute_c2v_from_c1_Picard and compute_log_act_coeffs_aq_chem already were):
!> - every concentration is paired with the valence / activity-coefficient parameters of ITS species
!>   through the water -> aqueous phase index maps (indices_aq_phase, ind_prim_species,
!>   ind_sec_species); the previous version read concentrations, valences and \f$\log\gamma_1\f$ by
!>   POSITION, which mixed species up whenever the aqueous species had been reordered (swapped
!>   primaries), e.g. an ionic strength of 1.47 instead of 0.66 for the saline water of the rezaei case
!> - the water activity coefficient (Raoult's law, compute_log_act_coeff_wat) is evaluated once after
!>   convergence instead of being refreshed inside the loop and fed back into the mass action law
!>
!> **Special Handling**:
!> - **Gas Phase**: the gas activity (fugacity) coefficients \f$\gamma_i = RT\f$ are constant and are
!>   computed once before the loop, for THIS gas object; the partial pressures \f$P_i=\gamma_i c_i\f$
!>   and the total pressure are evaluated after convergence
!> - **Non-ideal Solutions**: aqueous activity coefficients updated each iteration
!>
!> \note This version supports systems with aqueous + solid primary species (minerals, surface sites).
!>       For purely aqueous systems, use compute_c2_from_c1_aq_Picard instead.
!>
!> \warning Convergence is not guaranteed for all chemical systems. Poorly conditioned
!>          systems (extreme pH, high ionic strength) may require smaller tolerances or
!>          alternative methods (Newton-Raphson).
!>
!> \see compute_c2_from_c1_aq_Picard, compute_c2_from_c1_ideal, compute_c2v_from_c1_Picard
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2024
!>
!> \param[in,out] this Aqueous chemistry object containing chemical system and solution state
!> \param[in] c1 Primary species concentrations [M] (must be already set in the object)
!> \param[in] c2_ig Initial guess for secondary species concentrations [M]
!> \param[out] c2 Converged secondary species concentrations [M] (pre-allocated, size = num_sec_species)
!> \param[out] niter Number of Picard iterations performed until convergence or failure [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE if max iterations exceeded (e.g. the
!>             ionic strength of an absurd trial point makes the activity coefficients overflow): the
!>             object is then left untouched and c2 holds the last (non-converged, possibly NaN)
!>             iterate; the callers decide how to react (reject the point, fall back, or stop)
!>
subroutine compute_c2_from_c1_Picard(this,c1,c2_ig,c2,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this !< Aqueous chemistry object (polymorphic)
    real(kind=8), intent(in) :: c1(:) !< Primary concentrations [M] (dim=n_p)
    real(kind=8), intent(in) :: c2_ig(:) !< Initial guess secondary concentrations [M] (dim=n_e)
    real(kind=8), intent(out) :: c2(:) !< Secondary concentrations [M] (pre-allocated, dim=n_e)
    integer(kind=4), intent(out) :: niter !< Number of iterations [-]
    logical, intent(out) :: CV_flag !< TRUE if converges, FALSE otherwise
!> Local variables
    integer(kind=4) :: n_p,n_p_aq,n_2_aq,n_e,n_aq,n_diss,niter_max,i,j,k,m
    real(kind=8) :: rel_tol,log_rel_tol_lin,rel_change,log_rel_change,A_DH,I_fixed,ionic_str,lg_j
    real(kind=8), allocatable :: log_gamma1(:),log_gamma2(:),log_c2_cst(:),log_c2_new(:),log_c2_old(:),c2_old(:)
    real(kind=8), allocatable :: lg(:) !> log10 activity coefficients of the aqueous species at the current ionic strength, concentration order (n_aq)
    real(kind=8), allocatable :: coef_A(:),coef_beta(:),coef_gamma(:) !> per dissolved solid: -alpha*A*z^2, beta, gamma (n_diss)
    integer(kind=4), allocatable :: ind_diss(:) !> concentration index of each dissolved solid (n_diss)
    integer(kind=4), allocatable :: ind_prim_aq(:) !> concentration index of each aqueous primary species (n_p_aq)
    integer(kind=4), allocatable :: ind_sec_aq(:) !> concentration index of each aqueous secondary species, c2 order (n_2_aq)
    real(kind=8), allocatable :: z2_sec(:) !> squared valence of each aqueous secondary species (n_2_aq)
    logical, allocatable :: is_sec(:) !> concentration index belongs to an aqueous secondary species (n_aq)
!> Pre-process
    associate(sa => this%solid_chemistry%reactive_zone%speciation_alg, Se => this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star, &
              aq => this%aq_phase)
    n_p=sa%num_prim_species !< Total primary species
    n_p_aq=sa%num_aq_prim_species !< Aqueous primary species
    n_e=sa%num_eq_reactions !< Equilibrium reactions (= secondary species)
    n_2_aq=sa%num_sec_aq_species !< Secondary aqueous species
    n_aq=aq%num_species !< Total aqueous species
    n_diss=n_aq-aq%wat_flag !< Dissolved solids (water excluded)
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    !> log_rel_tol is stored as a base-10 logarithm by convention (e.g. -14 = 1e-14), matching the
    !> *_CV_params.dat files and the abs_tol/log_abs_tol pairing: convert it to a linear threshold
    log_rel_tol_lin=1d1**this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max
    !> Gas activity (fugacity) coefficients gamma_i = R*T are constant throughout the iteration, so
    !> they are computed once here, before get_log_gamma2 reads them. They must be set for THIS gas
    !> object (incl. per-cell target gases) because compute_partial_pressures reads them after convergence.
    if (associated(this%gas_chemistry)) call this%gas_chemistry%compute_log_act_coeffs_gases()
    !> Constant part of the mass action law (fixed along the iteration):
    !>   Se_1_star*log10(c1) + logK_tilde
    !>   + Se_1_star(:,non-aqueous primaries)*log gamma1(non-aqueous primaries)
    !>   - log gamma2(non-aqueous secondaries)
    log_gamma1=this%get_log_gamma1()
    log_gamma2=this%get_log_gamma2()
    allocate(log_c2_cst(n_e),log_c2_new(n_e),log_c2_old(n_e),c2_old(n_e))
    log_c2_cst=matmul(Se,log10(c1))+this%logK_tilde
    if (n_p>n_p_aq) log_c2_cst=log_c2_cst+matmul(Se(:,n_p_aq+1:n_p),log_gamma1(n_p_aq+1:n_p))
    if (n_e>n_2_aq) log_c2_cst(n_2_aq+1:n_e)=log_c2_cst(n_2_aq+1:n_e)-log_gamma2(n_2_aq+1:n_e)
    !> Concentration indices of the aqueous primary and secondary species (the latter in c2 order, as
    !> set_conc_sec_species scatters them) and the squared valences of the secondaries, the only ionic
    !> strength contributions that change along the iteration
    allocate(ind_prim_aq(n_p_aq),ind_sec_aq(n_2_aq),z2_sec(n_2_aq),is_sec(n_aq),lg(n_aq))
    ind_prim_aq=this%ind_prim_species(1:n_p_aq)
    ind_sec_aq=this%ind_sec_species(1:n_2_aq)
    is_sec=.false.
    do i=1,n_2_aq
        m=ind_sec_aq(i)
        is_sec(m)=.true.
        z2_sec(i)=dble(aq%aq_species(this%indices_aq_phase(m))%valence**2)
    end do
    !> Fixed part of the ionic strength: the aqueous species that are not secondary (primaries; each
    !> concentration paired with the valence of ITS species through the water -> aqueous phase index map)
    I_fixed=0d0
    do m=1,n_aq
        if (is_sec(m)) cycle
        I_fixed=I_fixed+this%concentrations(m)*dble(aq%aq_species(this%indices_aq_phase(m))%valence**2)
    end do
    I_fixed=0.5d0*I_fixed
    !> Activity coefficient model parameters of the dissolved solids, gathered once (water excluded:
    !> its log gamma is kept as stored in the object, like compute_log_act_coeffs_aq_chem does).
    !> beta and gamma are evaluated with THIS water's Debye-Hückel A and B (temperature dependent).
    allocate(ind_diss(n_diss),coef_A(n_diss),coef_beta(n_diss),coef_gamma(n_diss))
    A_DH=this%params_aq_sol%A
    do i=1,n_diss
        k=aq%ind_diss_solids(i)
        ind_diss(i)=this%indices_aq_species(k)
        coef_A(i)=-(aq%aq_species(k)%params_act_coeff%alpha*(A_DH*dble(aq%aq_species(k)%valence**2)))
        coef_beta(i)=aq%aq_species(k)%params_act_coeff%get_beta(this%params_aq_sol)
        coef_gamma(i)=aq%aq_species(k)%params_act_coeff%get_gamma(this%params_aq_sol)
    end do
    lg=this%log_act_coeffs !> non-dissolved-solid entries (water) stay as stored
!> Process: Picard iteration on the aqueous activity coefficients
    c2_old=c2_ig !< Initialize with user-provided initial guess
    log_c2_old=log10(c2_old) !< Log of initial guess for convergence check
    ionic_str=I_fixed+0.5d0*dot_product(z2_sec,c2_old(1:n_2_aq)) !> ionic strength of the initial guess
    call log_act_coeffs_diss_solids(ionic_str)
    niter=0 !< Initialize iteration counter to zero
    CV_flag=.false. !< Initialize convergence flag to FALSE (not converged yet)
    do
        niter=niter+1 !< Increment iteration counter
    !> Mass action law: constant part + aqueous activity coefficient corrections (aqueous primary
    !> columns of Se_1_star, then the aqueous secondaries' own gamma)
        log_c2_new=log_c2_cst
        do j=1,n_p_aq
            lg_j=lg(ind_prim_aq(j))
            if (lg_j/=0d0) log_c2_new=log_c2_new+Se(:,j)*lg_j
        end do
        do i=1,n_2_aq
            log_c2_new(i)=log_c2_new(i)-lg(ind_sec_aq(i))
        end do
        c2=1d1**log_c2_new !< Convert from log space to linear: c2 = 10^(log c2)
    !> Convergence: relative change of the concentrations AND of their logarithms (the latter against
    !> the linear threshold of log_rel_tol). A NaN iterate (overflow) is never accepted as converged.
        rel_change=0d0
        log_rel_change=0d0
        do i=1,n_e
            if (c2(i)/=c2(i)) then
                rel_change=huge(1d0)
                exit
            end if
            rel_change=max(rel_change,abs(c2(i)-c2_old(i))/abs(c2_old(i)))
            log_rel_change=max(log_rel_change,abs(log_c2_new(i)-log_c2_old(i))/abs(log_c2_old(i)))
        end do
    !> Ionic strength and aqueous activity coefficients at the new iterate: next iteration's input,
    !> or the converged state stored in the object below
        ionic_str=I_fixed+0.5d0*dot_product(z2_sec,c2(1:n_2_aq))
        call log_act_coeffs_diss_solids(ionic_str)
        if (rel_change<rel_tol .and. log_rel_change<log_rel_tol_lin) then
            CV_flag=.true. !< Converged successfully
            exit
        else if (niter==niter_max) then
            exit !> no CV: leave with CV_flag = .false. (the caller decides; the object is not updated)
        end if
        c2_old=c2 !< Update for next iteration convergence check
        log_c2_old=log_c2_new
    end do
    end associate
!> Post-process: state consistent with the converged c2 (skipped when the fixed point failed, so that
!> the object keeps the state of the last converged speciation; c2 then holds the last, possibly
!> NaN, iterate)
    if (CV_flag) then
    this%ionic_strength=ionic_str
    do i=1,n_diss
        this%log_act_coeffs(ind_diss(i))=lg(ind_diss(i)) !> log_10 activity coefficients of the dissolved solids at the converged ionic strength
    end do
    call this%set_conc_sec_species(c2) !> secondary concentrations in the aqueous, solid and gas chemistry objects
    call this%compute_activities() !< Final a_i = gamma_i * c_i for aqueous and solid phases
    call this%compute_pH()
    call this%compute_log_act_coeff_wat() !< log gamma_w via Raoult's law
    call this%compute_salinity() !< S = TDS / (1 + TDS) where TDS = sum(c_i * MW_i)
    call this%compute_alkalinity() !< Final alkalinity calculation
!> Final gas phase calculations (if present)
    if (associated(this%gas_chemistry)) then
        !> Partial pressures directly from the mass-action gas concentrations, P_i = gamma_i*c_i with
        !> gamma_i = R*T (computed once before the loop). Since logK_tilde uses the UNMODIFIED
        !> equilibrium constant (no partial-pressure fold), the mass-action c_i is the physical
        !> molarity (n/V) and P_i the correct equilibrium partial pressure.
        call this%gas_chemistry%compute_partial_pressures() !< P_i = gamma_i * c_i
        call this%gas_chemistry%compute_pressure() !< P_total = sum P_i
    end if
    end if
    deallocate(log_gamma1,log_gamma2,log_c2_cst,log_c2_new,log_c2_old,c2_old,lg,coef_A,coef_beta,coef_gamma,ind_diss, &
        ind_prim_aq,ind_sec_aq,z2_sec,is_sec)
contains
    !> log10 activity coefficients of the dissolved solids at ionic strength I (stored in lg by
    !> concentration index): -alpha*A*z^2*sqrt(I)/(1+beta*sqrt(I)) + gamma*I, as compute_log_act_coeffs_aq_phase
    subroutine log_act_coeffs_diss_solids(I)
        real(kind=8), intent(in) :: I
        real(kind=8) :: sqrt_I
        integer(kind=4) :: ii
        sqrt_I=sqrt(I)
        do ii=1,n_diss
            lg(ind_diss(ii))=coef_A(ii)*sqrt_I/(1d0+coef_beta(ii)*sqrt_I)+coef_gamma(ii)*I
        end do
    end subroutine log_act_coeffs_diss_solids
end subroutine compute_c2_from_c1_Picard
