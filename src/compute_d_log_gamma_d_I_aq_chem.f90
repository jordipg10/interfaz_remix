!> \file compute_d_log_gamma_d_I_aq_chem.f90
!> \brief Computes derivative of log₁₀(activity coefficients) with respect to ionic strength for aqueous variable activity species
!> \details This subroutine calculates the derivative ∂(log₁₀γᵢ)/∂I for all aqueous variable activity species 
!> (primary and secondary) using the extended Debye-Hückel activity coefficient model. This derivative is 
!> essential for computing Jacobian matrices in Newton-Raphson iterations with non-ideal activity coefficients.
!>
!> The extended Debye-Hückel equation for activity coefficients is:
!> \f[
!> \log_{10}(\gamma_i) = -\frac{\alpha_i A z_i^2 \sqrt{I}}{1 + \beta_i \sqrt{I}} + \gamma_i I
!> \f]
!> where:
!> - \f$ \gamma_i \f$ = activity coefficient of species i [-]
!> - \f$ \alpha_i \f$ = ion size parameter for species i [Å] or dimensionless scaling factor
!> - \f$ A \f$ = Debye-Hückel constant (temperature-dependent) [mol⁻¹/² kg¹/²]
!> - \f$ z_i \f$ = charge (valence) of species i [-]
!> - \f$ I \f$ = ionic strength [M]
!> - \f$ \beta_i \f$ = empirical parameter for species i [M⁻¹/²]
!> - \f$ \gamma_i \f$ = linear activity coefficient parameter (ion-ion interaction) [M⁻¹]
!>
!> The derivative is:
!> \f[
!> \frac{\partial \log_{10}(\gamma_i)}{\partial I} = -\frac{\alpha_i A z_i^2}{2\sqrt{I}(1 + \beta_i\sqrt{I})^2} + \gamma_i
!> \f]
!>
!> This derivative is used in computing activity coefficient corrections in the Jacobian for 
!> Newton-Raphson speciation algorithms.
!>
!> \param[in,out] this Aqueous chemistry object containing reactive zone, aqueous phase, ionic strength, and activity coefficient parameters
!> \param[inout] d_log_gamma_d_I Derivative of log₁₀(activity coefficients) w.r.t. ionic strength for all variable activity species (must be already allocated) [M⁻¹]

function compute_d_log_gamma_d_I_aq_chem(this) result(d_log_gamma_d_I)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    implicit none !> Enforce explicit variable declarations
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone, aqueous phase, ionic strength, and Debye-Hückel parameters [-]
    real(kind=8), allocatable :: d_log_gamma_d_I(:) !> Derivative of log₁₀(activity coefficients) w.r.t. ionic strength: ∂(log γ)/∂I (must be already allocated) [M⁻¹]

    integer(kind=4) :: i !> Loop counter for iterating over aqueous species [-]
    
    allocate(d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !> Allocate output array for all variable activity species (primary + secondary) [-]
    
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> Loop over all aqueous primary species
        d_log_gamma_d_I(i)=-this%aq_phase%aq_species(i)%params_act_coeff%alpha*(& !> Begin derivative computation: multiply by -αᵢ and continue on next line
            this%params_aq_sol%A*this%aq_phase%aq_species(i)%valence**2)/(2d0*sqrt(this%ionic_strength)*(& !> First term numerator: A·zᵢ², divided by 2√I times denominator part 1
            1d0+this%aq_phase%aq_species(i)%params_act_coeff%beta*sqrt(this%ionic_strength))**2) + & !> Denominator part 2: (1 + βᵢ√I)² for extended Debye-Hückel, then add linear term on next line
            this%aq_phase%aq_species(i)%params_act_coeff%gamma !> Add linear activity coefficient parameter γᵢ (ion-ion interaction term)
    end do !> End loop over aqueous primary species
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !> Loop over all aqueous secondary variable activity species
        d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=& !> Store derivative at index: num_primary + i (after all primary species in array)
            -this%aq_phase%aq_species(& !> Begin derivative computation: multiply by -αᵢ (accessing secondary species in aq_species array)
            this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%alpha*& !> Extract αᵢ for secondary aqueous species at index: num_aq_primary + i
            (this%params_aq_sol%A*this%aq_phase%aq_species(& !> First term numerator part 1: Debye-Hückel constant A
            this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%valence**2)/(2d0*sqrt(this%ionic_strength)& !> Numerator part 2: zᵢ², divided by 2√I
            *(1d0+this%aq_phase%aq_species(& !> Denominator part 1: (1 + βᵢ√I)², starting with 1 +
            this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%beta*& !> Extract βᵢ for secondary aqueous species
            sqrt(this%ionic_strength))**2) + & !> Complete denominator: βᵢ√I)², then add linear term on next line
            this%aq_phase%aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%gamma !> Add linear activity coefficient parameter γᵢ for secondary species
    end do !> End loop over aqueous secondary variable activity species
end function !> End of compute_d_log_gamma_d_I_aq_chem function