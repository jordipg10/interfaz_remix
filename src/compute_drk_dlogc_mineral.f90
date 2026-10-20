!> \file compute_drk_dlogc_mineral.f90
!> \brief Computes the gradient of a mineral kinetic reaction rate with respect to the LOGARITHM of the aqueous
!>        species concentrations (routine compute_drk_dlogc_mineral)
!> \details This subroutine calculates ∂rₖ/∂log₁₀c for a mineral dissolution/precipitation kinetic reaction.
!> The gradient is needed for constructing the Jacobian of the Newton residual of the reactive mixing
!> iterations, which is assembled directly in log-space (Newton_EI_eq_kin_anal(_ideal)_opt2 iterate on
!> p = log₁₀c₁): working with ∂/∂log₁₀c avoids forming 1/c_l factors that span as many orders of
!> magnitude as the concentrations themselves (and overflow for vanishing species), only to multiply them
!> back by c_l in the solver.
!>
!> The mineral kinetic reaction rate is given by a generalized rate law (see compute_rk_mineral):
!> \f[
!> r_k = \zeta A_s \exp\left(-\frac{E_a}{RT}\right) \sum_{j=1}^{n_{par}} k_j \prod_{i=1}^{n_{cat}} a_i^{p_{ji}} 
!> \left|\Omega_e^{\theta_j} - 1\right|^{\eta_j}, \qquad \zeta = \mathrm{sign}(\Omega_e-1),
!> \qquad \Omega_e = \min(\Omega,\Omega_{max})
!> \f]
!> where:
!> - \f$ A_s \f$ = reactive surface area [L^2/L^3]
!> - \f$ E_a \f$ = activation energy [J/mol]
!> - \f$ R \f$ = gas constant [J/(mol·K)]
!> - \f$ T \f$ = temperature [K]
!> - \f$ k_j \f$ = rate constant for j-th parallel reaction [mol/(L^2·T)]
!> - \f$ a_i \f$ = activity of i-th catalyser [-]
!> - \f$ p_{ji} \f$ = power for i-th catalyser in j-th reaction [-]
!> - \f$ \Omega \f$ = saturation index (Q/K) [-]
!> - \f$ \Omega_{max} \f$ = supersaturation threshold of the mineral (params%supersat_threshold; no clamp if <= 0) [-]
!> - \f$ \theta_j \f$ = reaction order parameter for j-th reaction [-]
!> - \f$ \eta_j \f$ = exponent parameter for j-th reaction [-]
!>
!> The gradient is computed assuming constant activity coefficients (\f$ \partial a/\partial \log_{10} c = \ln 10\, a \f$),
!> with \f$ \Omega = 10^{\log K}\prod_l a_l^{-\nu_l} \f$ (see compute_saturation_kin_min) and \f$ P_j = \prod_i a_i^{p_{ji}} \f$:
!> \f[
!> \frac{\partial r_k}{\partial \log_{10} c_l} = \ln 10\, A_s \exp\left(-\frac{E_a}{RT}\right) \left[
!> -\nu_l \sum_{j=1}^{n_{par}} k_j \eta_j \theta_j P_j \Omega_e^{\theta_j} \left|\Omega_e^{\theta_j} - 1\right|^{\eta_j-1}
!> + \zeta \sum_{j=1}^{n_{par}} k_j p_{jl} P_j \left|\Omega_e^{\theta_j} - 1\right|^{\eta_j} \right]
!> \f]
!> where \f$ \nu_l \f$ is the stoichiometric coefficient of species l and \f$ p_{jl} = 0 \f$ if l is not a catalyser.
!> The saturation term carries no \f$ \zeta \f$ because \f$ \frac{d}{d\Omega}\left[\zeta|\Omega^\theta-1|^\eta\right]
!> = \eta\theta\Omega^{\theta-1}|\Omega^\theta-1|^{\eta-1} \f$ (\f$ \zeta^2 = 1 \f$); the catalyser term does.
!> When the saturation index is clamped (\f$ \Omega > \Omega_{max} \f$) the rate no longer depends on Ω, so the
!> saturation term vanishes and only the catalyser term (evaluated at \f$ \Omega_e \f$) remains, consistently
!> with compute_rk_mineral. The arithmetic gradient, if needed, is \f$ \partial r_k/\partial c_l =
!> (\partial r_k/\partial \log_{10} c_l)/(\ln 10\, c_l) \f$.
!>
!> \param[in] this Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant
!> \param[in] conc_sp Concentrations of aqueous species participating in the reaction, ordered as indices_aq_phase [C]
!> \param[in] log_act_coeffs_sp Logarithm (base 10) of activity coefficients for aqueous species (gradient assumes gamma constant; used by the debug trace) [-]
!> \param[in] act_cat Activities of catalyser species (H⁺, OH⁻, etc.) [-]
!> \param[in] saturation Saturation index: Ω = Q/K (ion activity product over equilibrium constant), unclamped [-]
!> \param[in] react_surf Reactive surface area of mineral [m²]
!> \param[in] temp Temperature [K]
!> \param[out] drk_dlogc Gradient of reaction rate with respect to the log10 of the species concentrations: ∂rₖ/∂log₁₀c (must be already allocated) [mol/(kg·T)]
!>
!> \note [DRK-DBG] When the module switch dbg_kin_mineral (kin_mineral_m) is .true. — the solver arms it
!> only while the first target water is being solved during the first time step(s) — every input and
!> intermediate quantity (γ, a = γc, IAP rebuilt from them, Ω clamp, ζ, per-mechanism factors, both sums,
!> the unscaled and the final gradient) is written to the <root>_kin_mineral_dbg.out file (unit
!> dbg_kin_mineral_unit, opened by the solver; standard output if not open), one block per call, prefixed
!> with [DRK-DBG].

subroutine compute_drk_dlogc_mineral(this,conc_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dlogc)
    use kin_mineral_m, only: kin_mineral_c, dbg_kin_mineral, dbg_kin_mineral_unit, dbg_kin_mineral_step !> Import mineral kinetic reaction class and the [DRK-DBG] trace switch, unit and step
    use iso_fortran_env, only: output_unit !> Standard output: fallback of the [DRK-DBG] trace when its file is not open
    implicit none !> Enforce explicit variable declarations
    class(kin_mineral_c), intent(in) :: this !> Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant [-]
    real(kind=8), intent(in) :: conc_sp(:) !> Concentrations of aqueous species participating in the reaction, ordered as indices_aq_phase [C]
    real(kind=8), intent(in) :: log_act_coeffs_sp(:) !> Logarithm (base 10) of activity coefficients for aqueous species: log₁₀(γ) (gradient assumes γ constant; used by the debug trace) [-]
    real(kind=8), intent(in) :: act_cat(:) !> Activities of catalyser species (e.g., H⁺, OH⁻) that accelerate/inhibit reaction [-]
    real(kind=8), intent(in) :: saturation !> Saturation index: Ω = Q/K_eq (ratio of ion activity product to equilibrium constant), unclamped [-]
    real(kind=8), intent(in) :: react_surf !> Reactive surface area of mineral available for reaction [m²/kg water]
    real(kind=8), intent(in) :: temp !> Temperature (controls Arrhenius factor) [K]
    real(kind=8), intent(out) :: drk_dlogc(:) !> Gradient of reaction rate with respect to log10 species concentrations: ∂rₖ/∂log₁₀c (must be already allocated) [mol/(kg·T)]
    
    integer(kind=4) :: n_sp !> Number of aqueous species participating in the reaction [-]
    integer(kind=4) :: i !> Loop counter for catalyser species [-]
    integer(kind=4) :: l !> Loop counter for species with respect to which gradient is computed (current species) [-]
    integer(kind=4) :: j !> Loop counter for parallel reaction mechanisms [-]
    integer(kind=4) :: zeta !> Direction factor of compute_rk_mineral: -1 dissolution, +1 precipitation, 0 at Ω=1 [-]
    integer(kind=4), allocatable :: loc_cat(:) !> Position of each catalyser within conc_sp/drk_dlogc [-]
    real(kind=8), allocatable :: nu(:) !> Stoichiometric coefficient ν_l of each aqueous species of the reaction [-]
    real(kind=8) :: omega !> Saturation index actually used by the rate law: min(Ω, Ω_max) [-]
    logical :: clamped !> Ω > Ω_max: the rate is frozen at its Ω_max value, no dependence on Ω left
    real(kind=8) :: prod !> Product of catalyser activities: Pⱼ = ∏ᵢ aᵢ^pⱼᵢ [-]
    real(kind=8) :: sat_term !> |Ω^θⱼ - 1| [-]
    real(kind=8) :: term1 !> Contribution of mechanism j to sum1: kⱼ ηⱼ θⱼ Pⱼ Ω^θⱼ |Ω^θⱼ-1|^(ηⱼ-1) [mol/(L²·T)]
    real(kind=8) :: sum1 !> Saturation-term sum: Σⱼ kⱼ ηⱼ θⱼ Pⱼ Ω^θⱼ |Ω^θⱼ-1|^(ηⱼ-1) [mol/(L²·T)]
    real(kind=8), allocatable :: term2(:) !> Contribution of mechanism j to sum2 per catalyser: ζ kⱼ pⱼᵢ Pⱼ |Ω^θⱼ-1|^ηⱼ [mol/(L²·T)]
    real(kind=8), allocatable :: sum2(:) !> Catalyser-term sum per catalyser: ζ Σⱼ kⱼ pⱼᵢ Pⱼ |Ω^θⱼ-1|^ηⱼ [mol/(L²·T)]
    real(kind=8), allocatable :: grad_prod(:) !> ∂Pⱼ/∂log₁₀cᵢ / ln10 = pⱼᵢ·Pⱼ for each catalyser i [-]
    real(kind=8) :: arrh !> Arrhenius factor exp(-Eₐ/RT) [-]
    real(kind=8), parameter :: R=8.31446261815324 !> Universal gas constant [J/(mol·K)]
    real(kind=8), parameter :: ln10=log(1d1) !> ln(10): ∂a/∂log₁₀c = ln10·a [-]
    !> [DRK-DBG] trace variables (only evaluated when dbg_kin_mineral is .true.)
    integer(kind=4), save :: ncall_dbg=0 !> Number of traced calls so far, to tell successive Newton iterations/time steps apart [-]
    integer(kind=4) :: iu !> Output unit of the trace: the <root>_kin_mineral_dbg.out file if open, else standard output [-]
    real(kind=8) :: IAP_dbg !> ∏ₗ aₗ^(-νₗ) rebuilt from conc_sp and log_act_coeffs_sp, to check it against saturation [-]
    real(kind=8), allocatable :: gamma_dbg(:) !> γₗ = 10^log_act_coeffs_sp(l) [-]
    real(kind=8), allocatable :: act_dbg(:) !> aₗ = γₗ cₗ of the reaction's aqueous species [-]
    real(kind=8), allocatable :: ratio_dbg(:) !> act_cat(i)/a(loc_cat(i)): 1 if the caller's activities are consistent with its conc and log γ [-]
    
    iu=output_unit
    if (dbg_kin_mineral_unit/=-1) iu=dbg_kin_mineral_unit
    n_sp=size(conc_sp) !> Number of aqueous species of the reaction (= num_aq_rk)
    allocate(grad_prod(this%params%num_cat),sum2(this%params%num_cat),term2(this%params%num_cat),loc_cat(this%params%num_cat),&
        nu(n_sp))
    nu=this%stoichiometry(this%indices_react_species(1:n_sp))
    drk_dlogc=0d0
    sum1=0d0
    sum2=0d0
    !> cat_indices and indices_aq_phase are both aq_phase indices, so the catalyser slot in conc_sp/drk_dlogc
    !> is found by matching them (the caller's water-state indices are NOT valid here)
    do i=1,this%params%num_cat
        loc_cat(i)=0
        do l=1,n_sp
            if (this%indices_aq_phase(l)==this%params%cat_indices(i)) then
                loc_cat(i)=l
                exit
            end if
        end do
        if (loc_cat(i)==0) then
            write(0,*) 'compute_drk_dlogc_mineral: catalyser ', i, ' of reaction ', trim(this%name), &
                ' is not an aqueous species of the reaction; its derivative has no slot in drk_dlogc'
            error stop
        end if
    end do
    !> Supersaturation threshold (same clamp as compute_rk_mineral): beyond Ω_max the rate is frozen
    omega=saturation
    clamped=.false.
    if (this%params%supersat_threshold>0d0 .and. saturation>this%params%supersat_threshold) then
        omega=this%params%supersat_threshold
        clamped=.true.
    end if
    if (omega>1d0) then
        zeta=1
    else if (omega<1d0) then
        zeta=-1
    else
        zeta=0 !> catalyser term vanishes anyway (|Ω^θ-1|^η = 0)
    end if
    arrh=exp(-this%params%act_energy/(R*temp)) !> Arrhenius factor, as in compute_rk_mineral
    if (dbg_kin_mineral) then !> [DRK-DBG] inputs and the quantities derived from them
        ncall_dbg=ncall_dbg+1
        allocate(gamma_dbg(n_sp),act_dbg(n_sp),ratio_dbg(this%params%num_cat))
        gamma_dbg=1d1**log_act_coeffs_sp
        act_dbg=gamma_dbg*conc_sp
        ratio_dbg=-1d0 !> -1 = not computable (zero activity)
        where (act_dbg(loc_cat)>0d0) ratio_dbg=act_cat/act_dbg(loc_cat)
        if (all(act_dbg>0d0)) then
            IAP_dbg=product(act_dbg**(-nu))
        else
            IAP_dbg=-1d0 !> -1 = not computable (zero activity)
        end if
        write(iu,'(A)') '[DRK-DBG] ======================================================================'
        write(iu,'(A,A,A,I0,A,I0)') '[DRK-DBG] compute_drk_dlogc_mineral: reaction ', trim(this%name), ', time step ', &
            dbg_kin_mineral_step, ', traced call ', ncall_dbg
        write(iu,'(A,I0,A,I0,A,I0)') '[DRK-DBG]   n_sp = ', n_sp, ', num_par_reacts = ', this%params%num_par_reacts, &
            ', num_cat = ', this%params%num_cat
        write(iu,'(A,*(1X,I0))') '[DRK-DBG]   indices_aq_phase (aq. species of the reaction)   =', this%indices_aq_phase(1:n_sp)
        write(iu,'(A,*(1X,I0))') '[DRK-DBG]   indices_react_species                          =', this%indices_react_species(1:n_sp)
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   nu_l (stoichiometry of the aq. species)        =', nu
        write(iu,'(A,*(1X,I0))') '[DRK-DBG]   cat_indices (aq_phase index of each catalyser)  =', this%params%cat_indices
        write(iu,'(A,*(1X,I0))') '[DRK-DBG]   loc_cat (slot of each catalyser in conc_sp/drk_dlogc) =', loc_cat
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   conc_sp                        =', conc_sp
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   log_act_coeffs_sp              =', log_act_coeffs_sp
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   gamma_sp = 10**log_act_coeffs_sp =', gamma_dbg
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   act_sp = gamma_sp*conc_sp      =', act_dbg
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   act_cat (input)                =', act_cat
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   act_cat/act_sp(loc_cat) (1 if consistent, -1 n/a) =', ratio_dbg
        write(iu,'(A,ES14.6)') '[DRK-DBG]   saturation Omega (input)       = ', saturation
        write(iu,'(A,ES14.6)') '[DRK-DBG]   IAP = prod_l act_sp_l^(-nu_l)  = ', IAP_dbg
        if (IAP_dbg>0d0) write(iu,'(A,ES14.6)') '[DRK-DBG]   Omega/IAP (= 10**logK(T) used by the caller) = ', saturation/IAP_dbg
        write(iu,'(A,ES14.6,A,L1,A,ES14.6)') '[DRK-DBG]   Omega_max (supersat_threshold) = ', this%params%supersat_threshold, &
            ', clamped = ', clamped, ', Omega used = ', omega
        write(iu,'(A,I0)') '[DRK-DBG]   zeta = sign(Omega_used-1)      = ', zeta
        write(iu,'(A,ES14.6,A,ES14.6,A,ES14.6)') '[DRK-DBG]   react_surf = ', react_surf, ', temp = ', temp, &
            ', act_energy = ', this%params%act_energy
        write(iu,'(A,ES14.6,A,ES14.6)') '[DRK-DBG]   exp(-Ea/RT) = ', arrh, ', A_s*exp(-Ea/RT) = ', react_surf*arrh
        deallocate(gamma_dbg,act_dbg,ratio_dbg)
    end if
    do j=1,this%params%num_par_reacts
        prod=1d0
        do i=1,this%params%num_cat
            prod=prod*act_cat(i)**this%params%p(j,i) !> Pⱼ = ∏ᵢ aᵢ^pⱼᵢ
        end do
        grad_prod=this%params%p(j,:)*prod !> ∂Pⱼ/∂log₁₀cᵢ = ln10·pⱼᵢ·Pⱼ (γ constant); the ln10 is applied at the end
        sat_term=abs(omega**this%params%theta(j)-1d0)
        if (clamped) then
            term1=0d0 !> frozen rate: no dependence on Ω
        else
            term1=this%params%k(j)*prod*this%params%eta(j)*this%params%theta(j)*(omega**this%params%theta(j))*&
                sat_term**(this%params%eta(j)-1d0)
        end if
        term2=zeta*this%params%k(j)*grad_prod*sat_term**this%params%eta(j) !> rate carries ζ, so its catalyser derivative does too
        sum1=sum1+term1
        sum2=sum2+term2
        if (dbg_kin_mineral) then !> [DRK-DBG] factors of mechanism j
            write(iu,'(A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') '[DRK-DBG]   mechanism j = ', j, ': k_j = ', this%params%k(j), &
                ', theta_j = ', this%params%theta(j), ', eta_j = ', this%params%eta(j)
            write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]     p(j,:) (catalyser powers)                  =', this%params%p(j,:)
            write(iu,'(A,ES14.6)') '[DRK-DBG]     P_j = prod_i act_cat_i^p_ji                  = ', prod
            write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]     dP_j/dlog10c_i / ln10 = p_ji*P_j (per catalyser) =', grad_prod
            write(iu,'(A,ES14.6,A,ES14.6)') '[DRK-DBG]     Omega^theta_j = ', omega**this%params%theta(j), &
                ', |Omega^theta_j - 1| = ', sat_term
            write(iu,'(A,ES14.6)') '[DRK-DBG]     k_j*eta_j*theta_j*P_j*Omega^theta_j*|Omega^theta_j-1|^(eta_j-1) (-> sum1) = ', term1
            write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]     zeta*k_j*p_ji*P_j*|Omega^theta_j-1|^eta_j (-> sum2, per catalyser)  =', term2
        end if
    end do
    if (dbg_kin_mineral) then
        write(iu,'(A,ES14.6)') '[DRK-DBG]   sum1 (saturation-term sum)              = ', sum1
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   sum2 (catalyser-term sum, per catalyser) =', sum2
    end if
    do l=1,n_sp
        drk_dlogc(l)=-nu(l)*sum1 !> ∂Ω/∂log₁₀c_l = -ln10·ν_l·Ω (γ constant); the ln10 is applied at the end
    end do
    if (dbg_kin_mineral) write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   saturation term -nu_l*sum1               =', drk_dlogc
    do i=1,this%params%num_cat
        drk_dlogc(loc_cat(i))=drk_dlogc(loc_cat(i))+sum2(i)
    end do
    if (dbg_kin_mineral) write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   + catalyser term = drk_dlogc/(ln10 A_s exp(-Ea/RT)) =', drk_dlogc
    drk_dlogc=drk_dlogc*ln10*react_surf*arrh !> ln10·A_s·exp(-Eₐ/RT): A_s·exp(-Eₐ/RT) as in compute_rk_mineral, ln10 from ∂a/∂log₁₀c = ln10·a
    if (dbg_kin_mineral) then
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   drk_dlogc = dr/dlog10 c (output)         =', drk_dlogc
        write(iu,'(A,*(1X,ES14.6))') '[DRK-DBG]   (for reference dr/dc = drk_dlogc/(ln10 c) =', drk_dlogc/(ln10*conc_sp)
        write(iu,'(A)') '[DRK-DBG] ======================================================================'
    end if
    deallocate(grad_prod,sum2,term2,loc_cat,nu)
end subroutine