!> \file compute_drk_dc_mineral.f90
!> \brief Computes the gradient of mineral kinetic reaction rate with respect to aqueous species concentrations
!> \details This subroutine calculates ∂rₖ/∂c for a mineral dissolution/precipitation kinetic reaction.
!> The gradient is needed for constructing Jacobian matrices in reactive transport simulations with 
!> kinetic mineral reactions.
!>
!> The mineral kinetic reaction rate is given by a generalized rate law (see compute_rk_mineral):
!> \f[
!> r_k = \zeta A_s \exp\left(-\frac{E_a}{RT}\right) \sum_{j=1}^{n_{par}} k_j \prod_{i=1}^{n_{cat}} a_i^{p_{ji}} 
!> \left|\Omega^{\theta_j} - 1\right|^{\eta_j}, \qquad \zeta = \mathrm{sign}(\Omega-1)
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
!> - \f$ \theta_j \f$ = reaction order parameter for j-th reaction [-]
!> - \f$ \eta_j \f$ = exponent parameter for j-th reaction [-]
!>
!> The gradient is computed assuming constant activity coefficients (\f$ \partial a/\partial c = a/c \f$), with
!> \f$ \Omega = 10^{\log K}\prod_l a_l^{-\nu_l} \f$ (see compute_saturation_kin_min) and \f$ P_j = \prod_i a_i^{p_{ji}} \f$:
!> \f[
!> \frac{\partial r_k}{\partial c_l} = \frac{A_s \exp\left(-\frac{E_a}{RT}\right)}{c_l} \left[
!> -\nu_l \sum_{j=1}^{n_{par}} k_j \eta_j \theta_j P_j \Omega^{\theta_j} \left|\Omega^{\theta_j} - 1\right|^{\eta_j-1}
!> + \zeta \sum_{j=1}^{n_{par}} k_j p_{jl} P_j \left|\Omega^{\theta_j} - 1\right|^{\eta_j} \right]
!> \f]
!> where \f$ \nu_l \f$ is the stoichiometric coefficient of species l and \f$ p_{jl} = 0 \f$ if l is not a catalyser.
!> The saturation term carries no \f$ \zeta \f$ because \f$ \frac{d}{d\Omega}\left[\zeta|\Omega^\theta-1|^\eta\right]
!> = \eta\theta\Omega^{\theta-1}|\Omega^\theta-1|^{\eta-1} \f$ (\f$ \zeta^2 = 1 \f$); the catalyser term does.
!>
!> \param[in] this Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant
!> \param[in] conc_sp Concentrations of aqueous species participating in the reaction, ordered as indices_aq_phase [C]
!> \param[in] log_act_coeffs_sp Logarithm (base 10) of activity coefficients for aqueous species (unused: gamma assumed constant) [-]
!> \param[in] act_cat Activities of catalyser species (H⁺, OH⁻, etc.) [-]
!> \param[in] saturation Saturation index: Ω = Q/K (ion activity product over equilibrium constant) [-]
!> \param[in] react_surf Reactive surface area of mineral [m²]
!> \param[in] temp Temperature [K]
!> \param[out] drk_dc Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]

subroutine compute_drk_dc_mineral(this,conc_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dc)
    use kin_mineral_m, only: kin_mineral_c !> Import mineral kinetic reaction class
    implicit none !> Enforce explicit variable declarations
    class(kin_mineral_c), intent(in) :: this !> Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant [-]
    real(kind=8), intent(in) :: conc_sp(:) !> Concentrations of aqueous species participating in the reaction, ordered as indices_aq_phase [C]
    real(kind=8), intent(in) :: log_act_coeffs_sp(:) !> Logarithm (base 10) of activity coefficients for aqueous species: log₁₀(γ) (unused: γ assumed constant) [-]
    real(kind=8), intent(in) :: act_cat(:) !> Activities of catalyser species (e.g., H⁺, OH⁻) that accelerate/inhibit reaction [-]
    real(kind=8), intent(in) :: saturation !> Saturation index: Ω = Q/K_eq (ratio of ion activity product to equilibrium constant) [-]
    real(kind=8), intent(in) :: react_surf !> Reactive surface area of mineral available for reaction [m²/kg water]
    real(kind=8), intent(in) :: temp !> Temperature (controls Arrhenius factor) [K]
    real(kind=8), intent(out) :: drk_dc(:) !> Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]
    
    integer(kind=4) :: n_sp !> Number of aqueous species participating in the reaction [-]
    integer(kind=4) :: i !> Loop counter for catalyser species [-]
    integer(kind=4) :: l !> Loop counter for species with respect to which gradient is computed (current species) [-]
    integer(kind=4) :: j !> Loop counter for parallel reaction mechanisms [-]
    integer(kind=4) :: zeta !> Direction factor of compute_rk_mineral: -1 dissolution, +1 precipitation, 0 at Ω=1 [-]
    integer(kind=4), allocatable :: loc_cat(:) !> Position of each catalyser within conc_sp/drk_dc [-]
    real(kind=8) :: prod !> Product of catalyser activities: Pⱼ = ∏ᵢ aᵢ^pⱼᵢ [-]
    real(kind=8) :: sat_term !> |Ω^θⱼ - 1| [-]
    real(kind=8) :: sum1 !> Saturation-term sum: Σⱼ kⱼ ηⱼ θⱼ Pⱼ Ω^θⱼ |Ω^θⱼ-1|^(ηⱼ-1) [mol/(L²·T)]
    real(kind=8), allocatable :: sum2(:) !> Catalyser-term sum per catalyser: ζ Σⱼ kⱼ (pⱼᵢ Pⱼ/cᵢ) |Ω^θⱼ-1|^ηⱼ [mol/(L²·T·C)]
    real(kind=8), allocatable :: grad_prod(:) !> ∂Pⱼ/∂cᵢ = pⱼᵢ·Pⱼ/cᵢ for each catalyser i [1/C]
    real(kind=8), parameter :: R=8.31446261815324 !> Universal gas constant [J/(mol·K)]
    
    n_sp=size(conc_sp) !> Number of aqueous species of the reaction (= num_aq_rk)
    allocate(grad_prod(this%params%num_cat),sum2(this%params%num_cat),loc_cat(this%params%num_cat))
    drk_dc=0d0
    sum1=0d0
    sum2=0d0
    !> cat_indices and indices_aq_phase are both aq_phase indices, so the catalyser slot in conc_sp/drk_dc
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
            write(0,*) 'compute_drk_dc_mineral: catalyser ', i, ' of reaction ', trim(this%name), &
                ' is not an aqueous species of the reaction; its derivative has no slot in drk_dc'
            error stop
        end if
    end do
    if (saturation>1d0) then
        zeta=1
    else if (saturation<1d0) then
        zeta=-1
    else
        zeta=0 !> catalyser term vanishes anyway (|Ω^θ-1|^η = 0)
    end if
    do j=1,this%params%num_par_reacts
        prod=1d0
        do i=1,this%params%num_cat
            prod=prod*act_cat(i)**this%params%p(j,i) !> Pⱼ = ∏ᵢ aᵢ^pⱼᵢ
        end do
        do i=1,this%params%num_cat !> needs the complete Pⱼ, hence a separate loop
            grad_prod(i)=this%params%p(j,i)*prod/conc_sp(loc_cat(i)) !> ∂Pⱼ/∂cᵢ = pⱼᵢ·Pⱼ/cᵢ (γ constant)
        end do
        sat_term=abs(saturation**this%params%theta(j)-1d0)
        sum1=sum1+this%params%k(j)*prod*this%params%eta(j)*this%params%theta(j)*(saturation**this%params%theta(j))*&
            sat_term**(this%params%eta(j)-1d0)
        sum2=sum2+zeta*this%params%k(j)*grad_prod*sat_term**this%params%eta(j) !> rate carries ζ, so its catalyser derivative does too
    end do
    do l=1,n_sp
        drk_dc(l)=-this%stoichiometry(this%indices_react_species(l))*sum1/conc_sp(l) !> ∂Ω/∂c_l = -ν_l·Ω/c_l (γ constant)
    end do
    do i=1,this%params%num_cat
        drk_dc(loc_cat(i))=drk_dc(loc_cat(i))+sum2(i)
    end do
    drk_dc=drk_dc*react_surf*exp(-this%params%act_energy/(R*temp)) !> A_s·exp(-Eₐ/RT), as in compute_rk_mineral
    deallocate(grad_prod,sum2,loc_cat)
end subroutine