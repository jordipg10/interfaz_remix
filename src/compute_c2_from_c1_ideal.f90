!> \file compute_c2_from_c1_ideal.f90
!> \brief Computes secondary species concentrations from primary species concentrations using mass action law for ideal solutions
!> \details This subroutine explicitly computes all secondary species concentrations from primary species 
!> concentrations using the mass action law with ideal activity coefficient assumptions (γ = 1).
!> Unlike compute_c2_from_c1_aq_ideal, this routine accepts primary concentrations as input and can handle
!> both aqueous and solid primary species.
!>
!> Mass action law in logarithmic form:
!> \f[
!> \log_{10}(\mathbf{c}_2) = \mathbf{S}_{e,1}^* \cdot \log_{10}(\mathbf{c}_1) + \log_{10}\tilde{\mathbf{K}}
!> \f]
!> where:
!> - \f$ \mathbf{c}_2 \f$ = secondary species concentrations
!> - \f$ \mathbf{c}_1 \f$ = primary species concentrations (input parameter)
!> - \f$ \mathbf{S}_{e,1}^* \f$ = stoichiometric coefficient matrix (secondary vs. primary)
!> - \f$ \tilde{\mathbf{K}} \f$ = modified equilibrium constants
!>
!> Assumptions:
!> - Ideal activity coefficients: γᵢ = 1 for all species
!> - Equilibrium reactions are at equilibrium
!> - Primary concentrations provided as input (not from object)
!>
!> \param[in,out] this Aqueous chemistry object containing reactive zone and speciation algebra
!> \param[in] c1 Primary species concentrations: aqueous + solid primaries [C]
!> \param[out] c2 Secondary species concentrations (must be already allocated) [C]

subroutine compute_c2_from_c1_ideal(this,c1,log_act_coeffs,c2)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none !> Enforce explicit variable declarations
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone, aqueous phase, and speciation algebra [-]
    real(kind=8), intent(in) :: c1(:) !> Primary species concentrations: aqueous + solid (dimension = n_p) [C]
    real(kind=8), intent(in) :: log_act_coeffs(:) !> logarithm (base 10) of activity coefficients of species (0 for aqueous)
    real(kind=8), intent(out) :: c2(:) !> Secondary species concentrations (must be already allocated) [C]
    
    integer(kind=4) :: n_p !> Number of primary species (aqueous + solid) [-]
    integer(kind=4) :: n_sp !> Number of species
    real(kind=8), allocatable :: log_c2(:) !> Logarithm (base 10) of secondary species concentrations [-]
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Extract number of primary species from speciation algebra object
    n_sp=this%solid_chemistry%reactive_zone%speciation_alg%num_species !> Extract number of primary species from speciation algebra object
    
    log_c2=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,&
        log10(c1)+log_act_coeffs(1:n_p))+& !> Compute log₁₀(c₂) = Sₑ₁* · log₁₀(c₁) using stoichiometric coefficient matrix and logarithm of input primary concentrations
        this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde - & !> Add modified equilibrium constants log₁₀(K̃) to complete mass action law in logarithmic form
        log_act_coeffs(n_p+1:n_sp) !> Subtract log₁₀(γ₂) for secondary species from log₁₀(c₂)
    c2=10**log_c2 !> Convert from logarithmic form to linear concentrations: c₂ = 10^(log₁₀(c₂))
    call this%set_conc_sec_species(c2) !> Store all secondary species concentrations in aqueous chemistry object (both aqueous and other phases)
 end subroutine