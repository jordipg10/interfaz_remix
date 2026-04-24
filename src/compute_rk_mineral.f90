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
!>       \left| \Omega^{\theta_j} - 1 \right|^{\eta_j}
!> \f]
!> where:
!> - \f$A_s\f$ is the reactive surface area [m²]
!> - \f$E_a\f$ is the activation energy [J/mol]
!> - \f$R\f$ is the universal gas constant (8.31446 J/(mol·K))
!> - \f$T\f$ is the temperature [K]
!> - \f$\zeta\f$ is the direction: +1 for precipitation (Ω>1), -1 for dissolution (Ω<1) [-]
!> - \f$k_j\f$ is the rate constant for the j-th parallel mechanism [mol/(m²·s)]
!> - \f$a_i\f$ is the activity of the i-th catalyser [-]
!> - \f$p_{j,i}\f$ is the power/exponent for catalyser i in mechanism j [-]
!> - \f$\Omega\f$ is the saturation index (Q/K_eq) [-]
!> - \f$\theta_j\f$ is the reaction order for mechanism j [-]
!> - \f$\eta_j\f$ is the exponent parameter for mechanism j [-]
!>
!> \param[in] this Kinetic mineral object containing rate parameters (k, p, theta, eta, activation energy) [-]
!> \param[in] act_cat Activities of catalyser species (dimension = num_cat) [-]
!> \param[in] saturation Saturation index Ω = Q/K_eq (ion activity product / equilibrium constant) [-]
!> \param[in] react_surf Reactive surface area of the mineral [m²]
!> \param[in] temp Temperature [K]
!> \param[out] rk Kinetic reaction rate (positive for precipitation, negative for dissolution) [mol/s]

subroutine compute_rk_mineral(this,act_cat,saturation,react_surf,temp,rk)
    use kin_mineral_m, only: kin_mineral_c !> Import kinetic mineral class containing TST rate parameters
    implicit none !> Enforce explicit variable declarations
    class(kin_mineral_c), intent(in) :: this !> Kinetic mineral object containing rate parameters (k, p, theta, eta, num_par_reacts, num_cat, act_energy) [-]
    real(kind=8), intent(in) :: act_cat(:) !> Activities of catalyser species (dimension = num_cat, e.g., H⁺, OH⁻, ligands) [-]
    real(kind=8), intent(in) :: saturation !> Saturation index: Ω = Q/K_eq (ion activity product divided by equilibrium constant) [-]
    real(kind=8), intent(in) :: react_surf !> Reactive surface area of the mineral (specific surface area times volume fraction) [m²]
    real(kind=8), intent(in) :: temp !> Temperature (absolute temperature in Kelvin) [K]
    real(kind=8), intent(out) :: rk !> Kinetic reaction rate (output: positive for precipitation, negative for dissolution) [mol/s]
    
    integer(kind=4) :: zeta !> Direction indicator: +1 for precipitation (Ω>1), -1 for dissolution (Ω<1) [-]
    integer(kind=4) :: j !> Loop counter for parallel reaction mechanisms (j = 1 to num_par_reacts) [-]
    integer(kind=4) :: i !> Loop counter for catalyser species (i = 1 to num_cat) [-]
    real(kind=8) :: prod !> Product of catalyser activity terms: ∏(a_i^p_ji) for current mechanism j [-]
    real(kind=8), parameter :: R=8.31446261815324 !> Universal gas constant [J/(mol·K)]
    
    rk=0d0 !> Initialize kinetic reaction rate to zero before summing contributions from parallel mechanisms [mol/s]
    do j=1,this%params%num_par_reacts !> Loop over all parallel reaction mechanisms (e.g., H⁺-promoted, OH⁻-promoted, H₂O-promoted)
        prod=1d0 !> Initialize catalyser product term to 1 for mechanism j [-]
        do i=1,this%params%num_cat !> Loop over all catalyser species for this mechanism
            prod=prod*act_cat(i)**this%params%p(j,i) !> Multiply by catalyser activity raised to power p(j,i): prod *= a_i^p_ji [-]
        end do !> End catalyser loop
        rk=rk+prod*this%params%k(j)*abs((saturation**this%params%theta(j))-1d0)**this%params%eta(j) !> Add contribution of mechanism j: k_j · ∏(a_i^p_ji) · |Ω^θ_j - 1|^η_j [mol/s]
    end do !> End parallel mechanisms loop
    if (saturation>1d0) then !> Check if solution is supersaturated (Ω > 1, precipitation occurs)
        zeta=1 !> Set direction indicator to +1 for precipitation (mineral forms) [-]
    else if (saturation<1d0) then !> Check if solution is undersaturated (Ω < 1, dissolution occurs)
        zeta=-1 !> Set direction indicator to -1 for dissolution (mineral dissolves) [-]
    end if !> Note: if saturation = 1, zeta remains undefined (equilibrium, no net reaction)
    rk=rk*react_surf*exp(-this%params%act_energy/(R*temp)) !> Multiply by reactive surface area A_s and Arrhenius factor exp(-E_a/RT) [mol/s]
    rk=rk*zeta !> Apply direction indicator: positive rate for precipitation, negative for dissolution [mol/s]
end subroutine !> End of compute_rk_mineral subroutine