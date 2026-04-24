!> \file compute_drk_dc_Monod.f90
!> \brief Computes gradient of Monod kinetic reaction rate with respect to species concentrations
!> \details This subroutine calculates the gradient ∂rₖ/∂c for a redox kinetic reaction following 
!> Monod kinetics with inhibition, electron acceptor, and electron donor terms. The gradient is needed 
!> for constructing Jacobian matrices in reactive transport simulations with microbial reactions.
!>
!> The Monod kinetic rate law includes:
!> - Inhibition factors: \f$ \frac{K_{inh,j}}{K_{inh,j} + c_j} \f$
!> - Electron acceptor term: \f$ \frac{c_{acceptor}}{K_M + c_{acceptor}} \f$
!> - Electron donor term: \f$ \frac{c_{donor}}{K_M + c_{donor}} \f$
!>
!> The gradient has different forms depending on species type:
!> - For inhibitors:
!> \f[
!> \frac{\partial r_k}{\partial c_j} = -\frac{r_k}{K_{inh,j} + c_j}
!> \f]
!> - For electron acceptor/donor:
!> \f[
!> \frac{\partial r_k}{\partial c_M} = \frac{r_k \cdot K_M}{c_M(K_M + c_M)}
!> \f]
!>
!> Note: The implementation uses a simplified form with pre-computed rate rₖ.
!>
!> \param[in] this Redox kinetic reaction object containing Monod parameters (rate constant, half-saturation constants, inhibition constants)
!> \param[in] conc Species concentrations: conc = [conc_inh₁, ..., conc_inh_n, conc_acceptor, conc_donor] [C]
!> \param[in] rk Kinetic reaction rate (pre-computed) [M/T]
!> \param[out] drk_dc Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]

subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
    use redox_kin_reaction_m, only: redox_kin_c !> Import redox kinetic reaction class for Monod kinetics
    implicit none !> Enforce explicit variable declarations
    class(redox_kin_c), intent(in) :: this !> Redox kinetic reaction object containing Monod parameters (rate constant, K_M, K_inh) [-]
    real(kind=8), intent(in) :: conc(:) !> Species concentrations: conc = [conc_inh, conc_acceptor, conc_donor] where conc_inh may have multiple inhibitors [C]
    real(kind=8), intent(in) :: rk !> Kinetic reaction rate (pre-computed with all Monod terms) [M/T]
    real(kind=8), intent(out) :: drk_dc(:) !> Gradient of reaction rate with respect to species concentrations: ∂rₖ/∂c (must be already allocated) [T⁻¹]
    
    integer(kind=4) :: i !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: n !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: m !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: j !> Loop counter for inhibitors, electron acceptor, and electron donor [-]
    integer(kind=4) :: k !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: l !> Loop counter (not used in current implementation) [-]
    real(kind=8) :: prod_inh !> Product of inhibition factors (not used in current implementation) [-]
    real(kind=8), allocatable :: conc_inh(:) !> Inhibitor concentrations (not used in current implementation) [C]
    real(kind=8), allocatable :: conc_M(:) !> Electron acceptor and donor concentrations (not used in current implementation) [C]
    real(kind=8), allocatable :: prod_grad(:) !> Pre-factor for gradient computation: rₖ/c_M (temporary workaround for gradient calculation) [M/T/C]
    
    !prod_cat=1d0 !> Commented: Initialize product of catalyser terms to 1
    !prod_inh=1d0 !> Commented: Initialize product of inhibitor terms to 1
    allocate(prod_grad(2)) !> Allocate array for gradient pre-factors: one for electron acceptor, one for electron donor (temporary workaround)
    
    !rk=this%params%rate_cst !> Commented: Initialize reaction rate with rate constant k₀
    !prod_grad=this%params%rate_cst !> Commented: Initialize gradient pre-factor with rate constant
!> Inhibition factors (COMMENTED OUT - alternative approach to compute rate)
    !do j=1,this%params%n_inh !> Commented: Loop over all inhibitors
    !    rk=rk*this%params%k_inh(j)/(this%params%k_inh(j)+conc(j)) !> Commented: Multiply rate by inhibition factor: K_inh/(K_inh + c_inh)
    !end do
    !prod_grad=rk !> Commented: Store rate after inhibition for gradient calculation
    !rk=rk*prod_inh !> Commented: Multiply by product of all inhibition factors
!> Electron acceptor & donor (COMMENTED OUT - alternative approach)
    !do j=1,2 !> Commented: Loop over electron acceptor (j=1) and electron donor (j=2)
    !    rk=rk*conc(this%params%n_inh+j)/(this%params%k_M(j)+conc(this%params%n_inh+j)) !> Commented: Multiply by Monod term: c_M/(K_M + c_M)
    !end do
!> Temporary workaround for gradient pre-factor computation
    prod_grad(1)=prod_grad(1)*conc(this%params%n_inh+2)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+& !> Compute gradient pre-factor for electron acceptor: multiply by c_donor/(K_M,acceptor + c_acceptor)·(K_M,donor + c_donor)
        conc(this%params%n_inh+2))) !> Complete denominator for acceptor gradient pre-factor
    prod_grad(2)=prod_grad(2)*conc(this%params%n_inh+1)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+& !> Compute gradient pre-factor for electron donor: multiply by c_acceptor/(K_M,acceptor + c_acceptor)·(K_M,donor + c_donor)
        conc(this%params%n_inh+2))) !> Complete denominator for donor gradient pre-factor
    !rk=rk*prod_cat !> Commented: Multiply by product of catalyser terms (not currently used)
    !rk=rk*(1d0-conc(this%params%n_t+1)/this%params%cb_max) !> Commented: Multiply by logistic factor for biomass limitation: (1 - c_biomass/c_max)
!> Gradient computation: ∂rₖ/∂c for all species
    do j=1,this%params%n_inh !> Loop over all inhibitor species
        drk_dc(j)=-rk/(this%params%k_inh(j)+conc(j)) !> Compute gradient for j-th inhibitor: ∂rₖ/∂c_inh,j = -rₖ/(K_inh,j + c_j)
    end do !> End loop over inhibitors
    do j=1,2 !> Loop over electron acceptor (j=1) and electron donor (j=2)
        drk_dc(this%params%n_inh+j)=prod_grad(j)*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j)) !> Compute gradient for electron acceptor/donor: ∂rₖ/∂c_M = prod_grad·K_M/(K_M + c_M)
        !drk_dc(this%params%n_inh+j)=rk*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j))**2 !> Commented: Alternative formula: rₖ·K_M/(K_M + c_M)²
    end do !> End loop over electron acceptor and donor
end subroutine