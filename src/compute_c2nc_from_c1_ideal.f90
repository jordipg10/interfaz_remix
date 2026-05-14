!> \file compute_c2v_from_c1_ideal.f90
!> \brief Computes secondary variable activity species concentrations from primary species concentrations using mass action law for ideal solutions
!> \details This subroutine computes the concentrations of secondary variable activity species from the 
!> concentrations of primary species by explicitly applying the mass action law. The computation assumes 
!> ideal activity coefficients (γ = 1 for all species). This subroutine is designed for systems where 
!> primary species include both aqueous and solid species.
!>
!> For ideal solutions, activity equals concentration (a = c), so the mass action law becomes:
!> \f[
!> c_{2nc,k} = K_k^* \prod_{i=1}^{n_p} c_{1,i}^{S_{ki}^{nc*}}
!> \f]
!> where:
!> - \f$ c_{2nc,k} \f$ = concentration of secondary non-component variable activity species k
!> - \f$ K_k^* \f$ = modified equilibrium constant for species k
!> - \f$ c_{1,i} \f$ = concentration of primary species i
!> - \f$ S_{ki}^{nc*} \f$ = derived stoichiometric coefficient for non-component species
!> - \f$ n_p \f$ = number of primary species
!>
!> In logarithmic form:
!> \f[
!> \log_{10} c_{2nc,k} = \log_{10} K_k^* + \sum_{i=1}^{n_p} S_{ki}^{nc*} \log_{10} c_{1,i}
!> \f]
!> where \f$ S_{ki}^{nc*} = Se_{nc,1}^* \f$ is derived from stoichiometric matrix partitioning
!> for non-component secondary species.
!>
!> \param[in,out] this Aqueous chemistry object containing speciation algebra and equilibrium constants
!> \param[in] c1 Concentrations of primary species (aqueous + solid) [C]
!> \param[out] c2v Concentrations of secondary non-component variable activity species (must be already allocated) [C]

subroutine compute_c2v_from_c1_ideal(this,c1,log_act_coeffs,c2v)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none !> Enforce explicit variable declarations
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone and speciation algebra [-]
    real(kind=8), intent(in) :: c1(:) !> Concentrations of primary species (aqueous + solid, size = n_p) [C]
    real(kind=8), intent(in) :: log_act_coeffs(:) !> Logarithm (base 10) of activity coefficients of variable activity species (0 for aqueous species in ideal solution) [-]
    real(kind=8), intent(out) :: c2v(:) !> Concentrations of secondary non-component variable activity species (must be already allocated, size = n_e) [C]
    
    integer(kind=4) :: n_p,n_e !> n_p: number of primary species [-], n_v2_aq: number of aqueous secondary non-component variable activity species (unused) [-], n_e: number of equilibrium reactions [-]
    real(kind=8), allocatable :: log_c2v(:) !> log_c2v: log₁₀ of secondary non-component variable activity species concentrations [-]
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Extract number of primary species from speciation algebra [-]
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !> Extract number of equilibrium reactions from speciation algebra [-]

    if (n_p<0 .or. n_p>1000 .or. n_e<0 .or. n_e>1000) then
        print *, "DEBUG compute_c2v: CORRUPTED n_p=", n_p, " n_e=", n_e
        print *, "DEBUG compute_c2v: size(c1)=", size(c1), " size(log_act_coeffs)=", size(log_act_coeffs), " size(c2v)=", size(c2v)
        print *, "DEBUG compute_c2v: associated(solid_chemistry)=", associated(this%solid_chemistry)
        error stop "Corrupted speciation_alg in compute_c2v_from_c1_ideal"
    end if
    if (n_e==0) return !> No equilibrium reactions => nothing to compute
    !print *, this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star
    !print *, this%solid_chemistry%reactive_zone%speciation_alg%logK_star
    log_c2v=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,&
        log10(max(c1,1d-30))+log_act_coeffs(1:n_p))+& !> Compute log₁₀(c₂nc) = Se_nc_1_star·log₁₀(c₁) + logK_star (mass action law in logarithmic form for ideal solution where a = c) [-]
        this%solid_chemistry%reactive_zone%speciation_alg%logK_star - & !> Add log₁₀ of modified equilibrium constants (logK_star) for non-component secondary species [-]
        log_act_coeffs(n_p+1:) !> Subtract log₁₀ of activity coefficients for non-component secondary species (0 for ideal solution) [-]
    c2v=10**log_c2v !> Compute secondary non-component variable activity species concentrations by taking antilog: c₂nc = 10^(log₁₀(c₂nc)) [C]
    call this%set_conc_sec_var_act_species(c2v) !> Store computed secondary non-component variable activity species concentrations in the chemistry object [C]
 end subroutine !> End of compute_c2v_from_c1_ideal subroutine