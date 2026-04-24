!> \file compute_rk_drk_dc_Monod.f90
!> \brief Computes Monod kinetic reaction rate and its gradient with respect to species concentrations
!> \details This subroutine calculates the kinetic reaction rate for redox reactions using Monod kinetics,
!> which describes microbial growth and substrate utilization. The rate depends on:
!> - Inhibitor concentrations (reduce rate when present)
!> - Electron acceptor concentration (e.g., O₂, NO₃⁻, SO₄²⁻)
!> - Electron donor concentration (e.g., organic carbon, H₂)
!>
!> The Monod kinetic rate law is:
!> \f[
!> r_k = k_0 \cdot \prod_{j=1}^{n_{inh}} \frac{K_{inh,j}}{K_{inh,j} + c_{inh,j}} \cdot 
!>       \frac{c_{acc}}{K_{acc} + c_{acc}} \cdot \frac{c_{don}}{K_{don} + c_{don}}
!> \f]
!> where:
!> - \f$k_0\f$ is the maximum rate constant [M/(L³·T)]
!> - \f$K_{inh,j}\f$ is the inhibition half-saturation constant for inhibitor j [M/L³]
!> - \f$c_{inh,j}\f$ is the concentration of inhibitor j [M/L³]
!> - \f$K_{acc}\f$, \f$K_{don}\f$ are Monod half-saturation constants for acceptor and donor [M/L³]
!> - \f$c_{acc}\f$, \f$c_{don}\f$ are electron acceptor and donor concentrations [M/L³]
!>
!> The gradient is computed as:
!> - For inhibitors: \f$ \frac{\partial r_k}{\partial c_{inh}} = -\frac{r_k}{K_{inh} + c_{inh}} \f$
!> - For acceptor/donor: \f$ \frac{\partial r_k}{\partial c_M} = \frac{r_k \cdot K_M}{c_M (K_M + c_M)} \f$
!>
!> \param[in] this Redox kinetic reaction object containing Monod parameters (rate constant, half-saturation constants, inhibition constants) [-]
!> \param[in] conc Concentrations array: [inhibitors, electron acceptor, electron donor] [M/L³]
!> \param[out] rk Monod kinetic reaction rate [M/(L³·T)]
!> \param[out] drk_dc Gradient of reaction rate with respect to concentrations: ∂rk/∂c (must be pre-allocated) [T⁻¹]

subroutine compute_rk_drk_dc_Monod(this,conc,rk,drk_dc)
    use redox_kin_reaction_m, only: redox_kin_c !> Import redox kinetic reaction class containing Monod parameters
    implicit none !> Enforce explicit variable declarations
    class(redox_kin_c), intent(in) :: this !> Redox kinetic reaction object with Monod parameters (rate_cst, k_inh, k_M, n_inh) [-]
    real(kind=8), intent(in) :: conc(:) !> Concentrations array: conc = [c_inh(1:n_inh), c_acceptor, c_donor] [M/L³]
    real(kind=8), intent(out) :: rk !> Monod kinetic reaction rate (output) [M/(L³·T)]
    real(kind=8), intent(out) :: drk_dc(:) !> Gradient of reaction rate with respect to all concentrations: ∂rk/∂c (must be pre-allocated with size = n_inh + 2) [T⁻¹]
    
    integer(kind=4) :: i !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: n !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: m !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: j !> Loop counter for inhibitors and electron acceptor/donor [-]
    integer(kind=4) :: k !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: l !> Loop counter (not used in current implementation) [-]
    real(kind=8) :: prod_inh !> Product of inhibition factors (declared but not used in current implementation) [-]
    real(kind=8), allocatable :: conc_inh(:) !> Inhibitor concentrations (declared but not used in current implementation) [M/L³]
    real(kind=8), allocatable :: conc_M(:) !> Electron acceptor/donor concentrations (declared but not used in current implementation) [M/L³]
    real(kind=8), allocatable :: prod_grad(:) !> Pre-factor for gradient computation: rk before applying acceptor/donor terms (CHAPUZA - workaround for gradient calculation) [M/(L³·T)]
    
    !prod_cat=1d0 !> Commented out: product of catalyser terms (not used in Monod kinetics, would be for mineral kinetics) [-]
    !prod_inh=1d0 !> Commented out: product of inhibitor terms (computed directly in loop instead) [-]
    allocate(prod_grad(2)) !> Allocate gradient pre-factor array: size 2 for acceptor and donor (CHAPUZA - workaround for gradient computation) [-]
    
    rk=this%params%rate_cst !> Initialize reaction rate with maximum rate constant k₀ [M/(L³·T)]
    !prod_grad=this%params%rate_cst !> Commented out: alternative initialization of gradient pre-factor with rate constant
!> Inhibition factors: multiply rate by K_inh/(K_inh + c_inh) for each inhibitor
    do j=1,this%params%n_inh !> Loop over all inhibitors (n_inh = number of inhibitory species)
        rk=rk*this%params%k_inh(j)/(this%params%k_inh(j)+conc(j)) !> Apply inhibition factor: multiply rk by K_inh,j/(K_inh,j + c_inh,j) (reduces rate when inhibitor present) [M/(L³·T)]
    end do !> End inhibitor loop
    prod_grad=rk !> Store rate after inhibition factors (before acceptor/donor terms) for use in gradient calculation [M/(L³·T)]
    !rk=rk*prod_inh !> Commented out: alternative multiplication by pre-computed inhibition product
!> Electron acceptor & donor: apply Monod terms c_M/(K_M + c_M)
    do j=1,2 !> Loop over electron acceptor (j=1) and electron donor (j=2)
        rk=rk*conc(this%params%n_inh+j)/(this%params%k_M(j)+conc(this%params%n_inh+j)) !> Apply Monod term: multiply rk by c_M/(K_M + c_M) for acceptor or donor [M/(L³·T)]
    end do !> End acceptor/donor loop
!> CHAPUZA (workaround): Compute gradient pre-factors for acceptor and donor
    prod_grad(1)=prod_grad(1)*conc(this%params%n_inh+2)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+& !> Pre-factor for acceptor gradient: includes donor term c_don/(K_acc + c_acc)(K_don + c_don)
        conc(this%params%n_inh+2))) !> Denominator includes both acceptor and donor Monod terms [M/(L³·T)]
    prod_grad(2)=prod_grad(2)*conc(this%params%n_inh+1)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+& !> Pre-factor for donor gradient: includes acceptor term c_acc/(K_acc + c_acc)(K_don + c_don)
        conc(this%params%n_inh+2))) !> Denominator includes both acceptor and donor Monod terms [M/(L³·T)]
    !rk=rk*prod_cat !> Commented out: multiplication by catalyser product (not used in Monod kinetics)
    !rk=rk*(1d0-conc(this%params%n_t+1)/this%params%cb_max) !> Commented out: logistic growth factor for biomass limitation (1 - c_biomass/c_max)
!> Gradient computation
!> Gradient with respect to inhibitors: ∂rk/∂c_inh = -rk/(K_inh + c_inh)
    do j=1,this%params%n_inh !> Loop over all inhibitors
        drk_dc(j)=-rk/(this%params%k_inh(j)+conc(j)) !> Compute gradient for inhibitor j: negative because inhibitors decrease rate [T⁻¹]
    end do !> End inhibitor gradient loop
!> Gradient with respect to electron acceptor and donor: ∂rk/∂c_M = (rk·K_M)/(c_M·(K_M + c_M))
    do j=1,2 !> Loop over electron acceptor (j=1) and electron donor (j=2)
        drk_dc(this%params%n_inh+j)=prod_grad(j)*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j)) !> Compute gradient: (prod_grad · K_M)/(K_M + c_M), positive because acceptor/donor increase rate [T⁻¹]
        !drk_dc(this%params%n_inh+j)=rk*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j))**2 !> Commented out: alternative gradient formula rk·K_M/(K_M + c_M)² (mathematically equivalent after simplification)
    end do !> End acceptor/donor gradient loop
end subroutine !> End of compute_rk_drk_dc_Monod subroutine