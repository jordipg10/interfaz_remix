!> \file compute_rk_Monod.f90
!> \brief Computes Monod kinetic reaction rate for microbial redox reactions
!> \details This subroutine calculates the kinetic reaction rate for redox reactions using Monod kinetics,
!> which describes microbial growth and substrate utilization in biogeochemical systems. The rate law accounts for:
!> - Inhibitor concentrations (competitive or non-competitive inhibition)
!> - Electron acceptor concentration (terminal electron acceptor, e.g., O₂, NO₃⁻, SO₄²⁻)
!> - Electron donor concentration (substrate, e.g., organic carbon, H₂, CH₄)
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
!> - \f$K_{acc}\f$ is the Monod half-saturation constant for electron acceptor [M/L³]
!> - \f$K_{don}\f$ is the Monod half-saturation constant for electron donor [M/L³]
!> - \f$c_{acc}\f$ is the electron acceptor concentration [M/L³]
!> - \f$c_{don}\f$ is the electron donor concentration [M/L³]
!>
!> The Monod term c/(K+c) approaches 1 when c >> K (substrate abundant) and approaches c/K when c << K (substrate limiting).
!> The inhibition term K/(K+c) approaches 1 when c << K (no inhibition) and approaches 0 when c >> K (strong inhibition).
!>
!> \param[in] this Redox kinetic reaction object containing Monod parameters (rate_cst, k_inh, k_M, n_inh) [-]
!> \param[in] conc Concentrations array: [inhibitors, electron acceptor, electron donor] [M/L³]
!> \param[out] rk Monod kinetic reaction rate [M/(L³·T)]

subroutine compute_rk_Monod(this,conc,rk)
    use redox_kin_reaction_m, only: redox_kin_c !> Import redox kinetic reaction class containing Monod parameters
    implicit none !> Enforce explicit variable declarations
    class(redox_kin_c), intent(in) :: this !> Redox kinetic reaction object with Monod parameters (rate_cst, k_inh, k_M, n_inh) [-]
    real(kind=8), intent(in) :: conc(:) !> Concentrations array: conc = [c_inh(1:n_inh), c_acceptor, c_donor] [M/L³]
    real(kind=8), intent(out) :: rk !> Monod kinetic reaction rate (output) [M/(L³·T)]
    
    integer(kind=4) :: n !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: m !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: j !> Loop counter for inhibitors and electron acceptor/donor [-]
    integer(kind=4) :: k !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: l !> Loop counter (not used in current implementation) [-]
    real(kind=8) :: prod_cat !> Product of catalytic (substrate) terms for electron acceptor and donor [-]
    real(kind=8) :: prod_inh !> Product of inhibition terms for all inhibitors [-]
    real(kind=8), allocatable :: conc_inh(:) !> Inhibitor concentrations (declared but not used in current implementation) [M/L³]
    real(kind=8), allocatable :: conc_M(:) !> Electron acceptor/donor concentrations (declared but not used in current implementation) [M/L³]
    
    prod_cat=1d0 !> Initialize product of catalytic (substrate) terms to 1 (neutral multiplicative identity) [-]
    prod_inh=1d0 !> Initialize product of inhibition terms to 1 (no inhibition initially) [-]
    
    rk=this%params%rate_cst !> Initialize reaction rate with maximum rate constant k₀ [M/(L³·T)]
    do j=1,this%params%n_inh !> Loop over all inhibitors (n_inh = number of inhibitory species)
        prod_inh=prod_inh*this%params%k_inh(j)/(this%params%k_inh(j)+conc(j)) !> Accumulate inhibition factor: multiply by K_inh,j/(K_inh,j + c_inh,j) for inhibitor j [-]
    end do !> End inhibitor loop
    rk=rk*prod_inh !> Apply inhibition factors: multiply rate by product of all inhibition terms (reduces rate when inhibitors present) [M/(L³·T)]
    do j=1,2 !> Loop over electron acceptor (j=1) and electron donor (j=2)
        prod_cat=prod_cat*conc(this%params%n_inh+j)/(this%params%k_M(j)+conc(this%params%n_inh+j)) !> Accumulate Monod term: multiply by c_M/(K_M + c_M) for acceptor or donor [-]
    end do !> End acceptor/donor loop
    rk=rk*prod_cat !> Apply electron acceptor & donor Monod terms: multiply rate by product of substrate terms (increases rate with substrate availability) [M/(L³·T)]
    !rk=rk*(1d0-conc(this%params%n_t+1)/this%params%cb_max) !> Commented out: logistic growth factor for biomass limitation (1 - c_biomass/c_max), would limit rate when biomass approaches carrying capacity
end subroutine !> End of compute_rk_Monod subroutine