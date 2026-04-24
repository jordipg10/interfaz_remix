!> \file compute_drk_dc_mineral.f90
!> \brief Computes the gradient of mineral kinetic reaction rate with respect to aqueous species concentrations
!> \details This subroutine calculates ∂rₖ/∂c for a mineral dissolution/precipitation kinetic reaction.
!> The gradient is needed for constructing Jacobian matrices in reactive transport simulations with 
!> kinetic mineral reactions.
!>
!> The mineral kinetic reaction rate is given by a generalized rate law:
!> \f[
!> r_k = A_s \exp\left(-\frac{E_a}{RT}\right) \sum_{j=1}^{n_{par}} k_j \prod_{i=1}^{n_{cat}} a_i^{p_{ji}} 
!> \left[\left(\Omega^{\theta_j}\right)^{\eta_j} - 1\right]
!> \f]
!> where:
!> - \f$ A_s \f$ = reactive surface area [m²]
!> - \f$ E_a \f$ = activation energy [J/mol]
!> - \f$ R \f$ = gas constant [J/(mol·K)]
!> - \f$ T \f$ = temperature [K]
!> - \f$ k_j \f$ = rate constant for j-th parallel reaction [mol/(m²·s)]
!> - \f$ a_i \f$ = activity of i-th catalyser [-]
!> - \f$ p_{ji} \f$ = power for i-th catalyser in j-th reaction [-]
!> - \f$ \Omega \f$ = saturation index (Q/K) [-]
!> - \f$ \theta_j \f$ = reaction order parameter for j-th reaction [-]
!> - \f$ \eta_j \f$ = exponent parameter for j-th reaction [-]
!>
!> The gradient is computed as:
!> \f[
!> \frac{\partial r_k}{\partial c_l} = -\nu_l A_s \exp\left(-\frac{E_a}{RT}\right) \sum_{j=1}^{n_{par}} 
!> k_j \eta_j \theta_j \prod_i a_i^{p_{ji}} c_l^{\theta_j-1} \left(\prod_{m \neq l} a_m\right)^{\theta_j} 
!> \left(\Omega^{\theta_j} - 1\right)^{\eta_j-1}
!> \f]
!> where \f$ \nu_l \f$ is the stoichiometric coefficient of species l.
!>
!> \param[in] this Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant
!> \param[in] conc_sp Concentrations of aqueous species participating in the reaction [C]
!> \param[in] act_sp Activities of aqueous species participating in the reaction [-]
!> \param[in] log_act_coeffs_sp Logarithm (base 10) of activity coefficients for aqueous species [-]
!> \param[in] act_cat Activities of catalyser species (H⁺, OH⁻, etc.) [-]
!> \param[in] saturation Saturation index: Ω = Q/K (ion activity product over equilibrium constant) [-]
!> \param[in] react_surf Reactive surface area of mineral [m²]
!> \param[in] temp Temperature [K]
!> \param[out] drk_dc Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]

subroutine compute_drk_dc_mineral(this,conc_sp,act_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dc)
    use kin_mineral_m, only: kin_mineral_c !> Import mineral kinetic reaction class
    implicit none !> Enforce explicit variable declarations
    class(kin_mineral_c), intent(in) :: this !> Mineral kinetic reaction object containing reaction parameters, stoichiometry, and equilibrium constant [-]
    real(kind=8), intent(in) :: conc_sp(:) !> Concentrations of aqueous species participating in the reaction [C]
    real(kind=8), intent(in) :: act_sp(:) !> Activities of aqueous species participating in the reaction: a = γ·c [-]
    real(kind=8), intent(in) :: log_act_coeffs_sp(:) !> Logarithm (base 10) of activity coefficients for aqueous species: log₁₀(γ) [-]
    real(kind=8), intent(in) :: act_cat(:) !> Activities of catalyser species (e.g., H⁺, OH⁻) that accelerate/inhibit reaction [-]
    real(kind=8), intent(in) :: saturation !> Saturation index: Ω = Q/K_eq (ratio of ion activity product to equilibrium constant) [-]
    real(kind=8), intent(in) :: react_surf !> Reactive surface area of mineral available for reaction [m²]
    real(kind=8), intent(in) :: temp !> Temperature (controls Arrhenius factor) [K]
    real(kind=8), intent(out) :: drk_dc(:) !> Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]
    
    integer(kind=4) :: n_sp !> Number of aqueous species participating in the reaction [-]
    integer(kind=4) :: i !> Loop counter for catalyser species [-]
    integer(kind=4) :: l !> Loop counter for species with respect to which gradient is computed (current species) [-]
    integer(kind=4) :: j !> Loop counter for parallel reaction mechanisms [-]
    integer(kind=4) :: k !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: zeta !> Flag for dissolution (ζ=-1) or precipitation (ζ=1) (computed but not currently used) [-]
    real(kind=8) :: rk_j !> Reaction rate for j-th parallel mechanism (not used in current implementation) [M/T]
    real(kind=8) :: sum !> Sum over parallel reactions for gradient computation [T⁻¹]
    real(kind=8) :: prod !> Product of catalyser activities: ∏ᵢ aᵢ^pⱼᵢ [-]
    real(kind=8) :: prod_sat !> Product of activities excluding l-th species times activity coefficient and equilibrium constant: γₗ·K_eq·∏_{m≠l} aₘ [-]
    real(kind=8), parameter :: R=8.31446261815324 !> Universal gas constant [J/(mol·K)]
    
    drk_dc=0d0 !> Initialize gradient vector to zero before accumulation
    n_sp=size(act_sp) !> Extract number of aqueous species participating in reaction from size of activity array
    sum=0d0 !> Initialize sum over parallel reactions to zero
    !do l=1,n_sp !> Loop over all species to compute ∂rₖ/∂cₗ for each species l
    !    sum=0d0 !> Initialize sum over parallel reactions to zero for l-th species
    !    prod_sat=(10**log_act_coeffs_sp(l))*this%eq_cst !> Initialize saturation product: start with γₗ·K_eq (convert log₁₀(γₗ) to γₗ)
    !    do j=1,l-1 !> Loop over species indices before l-th species (j < l)
    !        prod_sat=prod_sat*act_sp(j) !> Multiply saturation product by activity of j-th species: accumulate ∏_{m<l} aₘ
    !    end do !> End loop over species before l
    !    do j=l+1,n_sp !> Loop over species indices after l-th species (j > l)
    !        prod_sat=prod_sat*act_sp(j) !> Multiply saturation product by activity of j-th species: accumulate ∏_{m>l} aₘ, giving γₗ·K_eq·∏_{m≠l} aₘ
    !    end do !> End loop over species after l
        do j=1,this%params%num_par_reacts !> Loop over all parallel reaction mechanisms (different rate laws for same mineral)
            prod=1d0 !> Initialize product of catalyser activities to 1 for j-th parallel reaction
            do i=1,this%params%num_cat !> Loop over all catalyser species (H⁺, OH⁻, etc.)
                prod=prod*act_cat(i)**this%params%p(j,i) !> Multiply by catalyser activity raised to power: accumulate ∏ᵢ aᵢ^pⱼᵢ
            end do !> End loop over catalysers
            sum=sum+prod*this%params%theta(j)*(saturation**this%params%theta(j))*(abs(&
                saturation**this%params%theta(j)-1d0)**(this%params%eta(j)-1d0))*&
                this%params%eta(j)*this%params%k(j)
            !sum=sum+prod*this%params%theta(j)*(conc_sp(l)**(this%params%theta(j)-1d0))*(prod_sat**this%params%theta(j))*& !> Add contribution from j-th parallel reaction: kⱼ·∏ᵢaᵢ^pⱼᵢ·θⱼ·cₗ^(θⱼ-1)·(γₗK∏_{m≠l}aₘ)^θⱼ, continue on next line
            !    this%params%k(j)*this%params%eta(j)*((& !> Multiply by rate constant kⱼ and exponent parameter ηⱼ, continue saturation term
            !    saturation**this%params%theta(j))-1d0)**(this%params%eta(j)-1d0) !> Multiply by (Ω^θⱼ - 1)^(ηⱼ-1) for TST-based rate law
        end do !> End loop over parallel reactions
        !drk_dc(l)=sum*(-this%stoichiometry(l)) !> Compute gradient component: multiply sum by negative stoichiometric coefficient -νₗ
    !end do !> End loop over species
    do l=1,n_sp !> Loop over all species to compute ∂rₖ/∂cₗ for each species l
        drk_dc(l)=sum*(-this%stoichiometry(this%indices_react_species(l)))/conc_sp(l) !> Compute gradient component: multiply sum by negative stoichiometric coefficient -νₗ
    end do !> End loop over species
!> Flag for dissolution or precipitation (computed but not currently used in gradient)
    !if (saturation<1d0) then !> Check if saturation index less than 1 (undersaturated condition)
    !    zeta=-1 !> Set flag to -1 for dissolution (mineral dissolves when undersaturated)
    !else if (saturation>1d0) then !> Check if saturation index greater than 1 (supersaturated condition)
    !    zeta=1 !> Set flag to +1 for precipitation (mineral precipitates when supersaturated)
    !end if !> End saturation check (zeta undefined when Ω=1, equilibrium)
    drk_dc=drk_dc*react_surf*exp(-this%params%act_energy/(R*temp)) !> Multiply gradient by reactive surface area and Arrhenius factor: exp(-Eₐ/RT)
    !> Apply dissolution/precipitation sign (consistent with compute_rk_mineral)
    ! if (saturation>1d0) then
    !     zeta=1
    ! else if (saturation<1d0) then
    !     zeta=-1
    ! else
    !     zeta=0
    ! end if
    ! drk_dc=drk_dc*zeta
end subroutine