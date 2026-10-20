!> \file compute_rk_mineral.f90
!> \brief Computes mineral kinetic reaction rate using transition state theory (TST) based rate law
!> \details This subroutine calculates the kinetic reaction rate for mineral dissolution or precipitation
!> using a generalized TST-based rate law that accounts for:
!> - Parallel reaction mechanisms (multiple pathways with different catalysers)
!> - Catalyser effects (e.g., H⁺, OH⁻, organic ligands)
!> - Saturation state dependence (thermodynamic drive)
!> - Reactive surface area
!> - Temperature dependence via Arrhenius equation
!>
!> The general rate law is:
!> \f[
!> r_k = A_s \cdot \exp\left(-\frac{E_a}{RT}\right) \cdot \zeta \cdot 
!>       \sum_{j=1}^{n_{par}} k_j \cdot \prod_{i=1}^{n_{cat}} a_i^{p_{j,i}} \cdot 
!>       \left| \Omega_e^{\theta_j} - 1 \right|^{\eta_j}, \qquad \Omega_e = \min(\Omega,\Omega_{max})
!> \f]
!> where:
!> - \f$A_s\f$ is the reactive surface area [m²]
!> - \f$E_a\f$ is the activation energy [J/mol]
!> - \f$R\f$ is the universal gas constant (8.31446 J/(mol·K))
!> - \f$T\f$ is the temperature [K]
!> - \f$\zeta\f$ is the direction: +1 for precipitation (Ω_e>1), -1 for dissolution (Ω_e<1), 0 at equilibrium [-]
!> - \f$k_j\f$ is the rate constant for the j-th parallel mechanism [mol/(m²·s)]
!> - \f$a_i\f$ is the activity of the i-th catalyser [-]
!> - \f$p_{j,i}\f$ is the power/exponent for catalyser i in mechanism j [-]
!> - \f$\Omega\f$ is the saturation index (Q/K_eq) [-]
!> - \f$\Omega_{max}\f$ is the supersaturation threshold of the mineral (params%supersat_threshold, read from the
!>   kinetics database): beyond it the precipitation rate is frozen at its Ω_max value (continuous clamp, zero
!>   slope in Ω, mirrored by compute_drk_dlogc_mineral). Ω_max <= 0 means no clamp; Ω_max = 1 means no precipitation [-]
!> - \f$\theta_j\f$ is the reaction order for mechanism j [-]
!> - \f$\eta_j\f$ is the exponent parameter for mechanism j [-]
!>
!> \param[in] this Kinetic mineral object containing rate parameters (k, p, theta, eta, activation energy) [-]
!> \param[in] act_cat Activities of catalyser species (dimension = num_cat) [-]
!> \param[in] saturation Saturation index Ω = Q/K_eq (ion activity product / equilibrium constant), unclamped [-]
!> \param[in] react_surf Reactive surface area of the mineral [m²]
!> \param[in] temp Temperature [K]
!> \param[out] rk Kinetic reaction rate (positive for precipitation, negative for dissolution) [mol/s]
!>
!> \note [RK-DBG] When the module switch dbg_kin_mineral (kin_mineral_m) is .true. — the solver arms it only
!> while the first target water is being solved during the first time step(s) — every input and intermediate
!> quantity (Ω clamp, ζ, Arrhenius factor, per-mechanism factors, the mechanism sum and the final rate) is
!> written to the <root>_kin_mineral_dbg.out file (unit dbg_kin_mineral_unit, opened by the solver; standard
!> output if not open), one block per call, prefixed with [RK-DBG].

subroutine compute_rk_mineral(this,act_cat,saturation,react_surf,temp,rk)
    use kin_mineral_m, only: kin_mineral_c, dbg_kin_mineral, dbg_kin_mineral_unit, dbg_kin_mineral_step !> Import kinetic mineral class containing TST rate parameters and the [RK-DBG] trace switch, unit and step
    use iso_fortran_env, only: output_unit !> Standard output: fallback of the [RK-DBG] trace when its file is not open
    implicit none !> Enforce explicit variable declarations
    class(kin_mineral_c), intent(in) :: this !> Kinetic mineral object containing rate parameters (k, p, theta, eta, num_par_reacts, num_cat, act_energy) [-]
    real(kind=8), intent(in) :: act_cat(:) !> Activities of catalyser species (dimension = num_cat, e.g., H⁺, OH⁻, ligands) [-]
    real(kind=8), intent(in) :: saturation !> Saturation index: Ω = Q/K_eq (ion activity product divided by equilibrium constant) [-]
    real(kind=8), intent(in) :: react_surf !> Reactive surface area of the mineral (specific surface area) [m²/kg_water]
    real(kind=8), intent(in) :: temp !> Temperature (absolute temperature in Kelvin) [K]
    real(kind=8), intent(out) :: rk !> Kinetic reaction rate (output: positive for precipitation, negative for dissolution)
    
    integer(kind=4) :: zeta !> Direction indicator: +1 for precipitation (Ω>1), -1 for dissolution (Ω<1), 0 at equilibrium (Ω=1) [-]
    integer(kind=4) :: j !> Loop counter for parallel reaction mechanisms (j = 1 to num_par_reacts) [-]
    integer(kind=4) :: i !> Loop counter for catalyser species (i = 1 to num_cat) [-]
    real(kind=8) :: omega !> Saturation index actually used by the rate law: min(Ω, Ω_max) [-]
    logical :: clamped !> Ω > Ω_max: the rate is frozen at its Ω_max value
    real(kind=8) :: prod !> Product of catalyser activity terms: ∏(a_i^p_ji) for current mechanism j [-]
    real(kind=8) :: sat_term !> |Ω^θⱼ - 1| for current mechanism j [-]
    real(kind=8) :: term !> Contribution of mechanism j: kⱼ·Pⱼ·|Ω^θⱼ-1|^ηⱼ [mol/(m²·s)]
    real(kind=8) :: arrh !> Arrhenius factor exp(-Eₐ/RT) [-]
    real(kind=8), parameter :: R=8.31446261815324 !> Universal gas constant [J/(mol·K)]
    integer(kind=4), save :: ncall_dbg=0 !> [RK-DBG] number of traced calls so far, to tell successive Newton iterations/time steps apart [-]
    integer(kind=4) :: iu !> [RK-DBG] output unit of the trace: the <root>_kin_mineral_dbg.out file if open, else standard output [-]
    
    iu=output_unit
    if (dbg_kin_mineral_unit/=-1) iu=dbg_kin_mineral_unit
    !> Supersaturation threshold (same clamp as compute_drk_dlogc_mineral): beyond Ω_max the rate is frozen
    omega=saturation
    clamped=.false.
    if (this%params%supersat_threshold>0d0 .and. saturation>this%params%supersat_threshold) then
        omega=this%params%supersat_threshold
        clamped=.true.
    end if
    if (omega>1d0) then !> Check if solution is supersaturated (Ω > 1, precipitation occurs)
        zeta=1 !> Set direction indicator to +1 for precipitation (mineral forms) [-]
    else if (omega<1d0) then !> Check if solution is undersaturated (Ω < 1, dissolution occurs)
        zeta=-1 !> Set direction indicator to -1 for dissolution (mineral dissolves) [-]
    else
        zeta=0 !> Equilibrium (Ω = 1): no net reaction (the mechanism sum is 0 anyway, |Ω^θ-1|^η = 0)
    end if
    arrh=exp(-this%params%act_energy/(R*temp)) !> Arrhenius factor exp(-E_a/RT) [-]
    if (dbg_kin_mineral) then !> [RK-DBG] inputs and the quantities derived from them
        ncall_dbg=ncall_dbg+1
        write(iu,'(A)') '[RK-DBG] ======================================================================'
        write(iu,'(A,A,A,I0,A,I0)') '[RK-DBG] compute_rk_mineral: reaction ', trim(this%name), ', time step ', dbg_kin_mineral_step, &
            ', traced call ', ncall_dbg
        write(iu,'(A,I0,A,I0)') '[RK-DBG]   num_par_reacts = ', this%params%num_par_reacts, ', num_cat = ', this%params%num_cat
        write(iu,'(A,*(1X,I0))') '[RK-DBG]   cat_indices (aq_phase index of each catalyser) =', this%params%cat_indices
        write(iu,'(A,*(1X,ES14.6))') '[RK-DBG]   act_cat (input)                =', act_cat
        write(iu,'(A,ES14.6)') '[RK-DBG]   saturation Omega (input)       = ', saturation
        write(iu,'(A,ES14.6,A,L1,A,ES14.6)') '[RK-DBG]   Omega_max (supersat_threshold) = ', this%params%supersat_threshold, &
            ', clamped = ', clamped, ', Omega used = ', omega
        write(iu,'(A,I0)') '[RK-DBG]   zeta = sign(Omega_used-1)      = ', zeta
        write(iu,'(A,ES14.6,A,ES14.6,A,ES14.6)') '[RK-DBG]   react_surf = ', react_surf, ', temp = ', temp, &
            ', act_energy = ', this%params%act_energy
        write(iu,'(A,ES14.6,A,ES14.6)') '[RK-DBG]   exp(-Ea/RT) = ', arrh, ', A_s*exp(-Ea/RT) = ', react_surf*arrh
    end if
    rk=0d0 !> Initialize kinetic reaction rate to zero before summing contributions from parallel mechanisms [mol/s]
    do j=1,this%params%num_par_reacts !> Loop over all parallel reaction mechanisms (e.g., H⁺-promoted, OH⁻-promoted, H₂O-promoted)
        prod=1d0 !> Initialize catalyser product term to 1 for mechanism j [-]
        do i=1,this%params%num_cat !> Loop over all catalyser species for this mechanism
            prod=prod*act_cat(i)**this%params%p(j,i) !> Multiply by catalyser activity raised to power p(j,i): prod *= a_i^p_ji [-]
        end do !> End catalyser loop
        sat_term=abs((omega**this%params%theta(j))-1d0) !> |Ω_e^θ_j - 1| [-]
        term=prod*this%params%k(j)*sat_term**this%params%eta(j) !> Contribution of mechanism j: k_j · ∏(a_i^p_ji) · |Ω_e^θ_j - 1|^η_j [mol/s]
        rk=rk+term !> Add contribution of mechanism j
        if (dbg_kin_mineral) then !> [RK-DBG] factors of mechanism j
            write(iu,'(A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') '[RK-DBG]   mechanism j = ', j, ': k_j = ', this%params%k(j), &
                ', theta_j = ', this%params%theta(j), ', eta_j = ', this%params%eta(j)
            write(iu,'(A,*(1X,ES14.6))') '[RK-DBG]     p(j,:) (catalyser powers)                  =', this%params%p(j,:)
            write(iu,'(A,ES14.6)') '[RK-DBG]     P_j = prod_i act_cat_i^p_ji                  = ', prod
            write(iu,'(A,ES14.6,A,ES14.6,A,ES14.6)') '[RK-DBG]     Omega^theta_j = ', omega**this%params%theta(j), &
                ', |Omega^theta_j - 1| = ', sat_term, ', |Omega^theta_j - 1|^eta_j = ', sat_term**this%params%eta(j)
            write(iu,'(A,ES14.6)') '[RK-DBG]     k_j*P_j*|Omega^theta_j-1|^eta_j (-> sum)       = ', term
        end if
    end do !> End parallel mechanisms loop
    if (dbg_kin_mineral) write(iu,'(A,ES14.6)') '[RK-DBG]   sum_j k_j*P_j*|Omega^theta_j-1|^eta_j     = ', rk
    rk=rk*react_surf*arrh !> Multiply by reactive surface area A_s and Arrhenius factor exp(-E_a/RT) [mol/s]
    rk=rk*zeta !> Apply direction indicator: positive rate for precipitation, negative for dissolution [mol/s]
    if (dbg_kin_mineral) then
        write(iu,'(A,ES14.6)') '[RK-DBG]   rk = zeta*A_s*exp(-Ea/RT)*sum (output)    = ', rk
        write(iu,'(A)') '[RK-DBG] ======================================================================'
    end if
end subroutine !> End of compute_rk_mineral subroutine