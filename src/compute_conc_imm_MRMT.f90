!> @file compute_conc_imm_MRMT.f90
!> @brief Computes concentrations in immobile zones for Multi-Rate Mass Transfer (MRMT) model
!> @details This subroutine updates concentrations in immobile zones at time step k+1 using a
!> theta-weighted time integration scheme. The MRMT model represents physical or chemical
!> heterogeneity by partitioning the porous medium into mobile and multiple immobile zones,
!> each with distinct exchange rates.
!>
!> @par Mathematical Formulation:
!> For each immobile zone j, the concentration update follows:
!> 
!> c_imm^(k+1) = [c_imm^k * (1 - α*Δt*(1-θ)) + α*Δt*(θ*c_m^(k+1) + (1-θ)*c_m^k)] / (1 + α*Δt*θ)
!>
!> where:
!>   - c_imm = concentration in immobile zone [M/L³]
!>   - c_m = concentration in mobile zone [M/L³]
!>   - α = first-order exchange rate coefficient [1/T]
!>   - Δt = time step size [T]
!>   - θ = time weighting factor (0 ≤ θ ≤ 1):
!>     * θ = 0: explicit scheme (forward Euler)
!>     * θ = 0.5: Crank-Nicolson (trapezoidal rule)
!>     * θ = 1: fully implicit scheme (backward Euler)
!>
!> @par Physical Interpretation:
!> The exchange rate α characterizes the rate of mass transfer between mobile and immobile zones.
!> Fast exchange (large α) approaches local equilibrium, while slow exchange (small α) creates
!> non-equilibrium behavior with tailing in breakthrough curves.
!>
!> @see MRMT_m Multi-Rate Mass Transfer module
!> @see Haggerty & Gorelick (1995), Water Resources Research, 31(10), 2383-2400
!> @author Generated documentation
!> @date November 2025

!> @brief Update immobile zone concentrations at next time step using theta-weighted scheme
!> @param[in] this MRMT object containing immobile zone properties (exchange rates, number of zones)
!> @param[in] theta Time weighting factor θ ∈ [0,1] for temporal discretization (0=explicit, 0.5=Crank-Nicolson, 1=implicit) [-]
!> @param[in] conc_imm_old Concentrations in immobile zones at time k: c_imm^k [M/L³]
!> @param[in] conc_mob_old Concentration in mobile zone at time k: c_m^k [M/L³]
!> @param[in] conc_mob_new Concentration in mobile zone at time k+1: c_m^(k+1) [M/L³]
!> @param[in] Delta_t Time step size Δt [T]
!> @param[out] conc_imm_new Concentrations in immobile zones at time k+1: c_imm^(k+1) [M/L³]
subroutine compute_conc_imm_MRMT(this,theta,conc_imm_old,conc_mob_old,conc_mob_new,Delta_t,conc_imm_new)
    use MRMT_m, only: MRMT_c !< Import Multi-Rate Mass Transfer class definition
    implicit none !< Enforce explicit variable declarations for type safety
    
    class(MRMT_c), intent(in) :: this !< MRMT object containing immobile zone parameters (n_imm zones, exchange rates α) [-]
    real(kind=8), intent(in) :: theta !< Time weighting factor θ: 0 (explicit) ≤ θ ≤ 1 (implicit) [-]
    real(kind=8), intent(in) :: conc_imm_old(:) !< Immobile zone concentrations at time k: c_imm^k [M/L³]
    real(kind=8), intent(in) :: conc_mob_old(:) !< Mobile zone concentration at time k: c_m^k [M/L³]
    real(kind=8), intent(in) :: conc_mob_new(:) !< Mobile zone concentration at time k+1: c_m^(k+1) [M/L³]
    real(kind=8), intent(in) :: Delta_t !< Time step size Δt [T]
    real(kind=8), intent(out) :: conc_imm_new(:) !< Output: immobile zone concentrations at time k+1: c_imm^(k+1) [M/L³]
    
    integer(kind=8) :: j !< Loop counter for iterating over immobile zones (1 to n_imm) [-]
    
    !> Validate time weighting factor is within physically meaningful range
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1" !< Ensure 0 ≤ θ ≤ 1 (stability and consistency requirement)
    
    !> Loop over all immobile zones to update each concentration independently
    do j=1,this%n_imm !< Iterate through all immobile zones (each has distinct exchange rate α_j)
        !> Update concentration in immobile zone j using theta-weighted finite difference scheme
        !> Numerator terms:
        !>   1. c_imm^k * (1 - α*Δt*(1-θ)): old concentration diminished by explicit exchange
        !>   2. α*Δt*(θ*c_m^(k+1) + (1-θ)*c_m^k): influx from mobile zone (theta-weighted)
        !> Denominator: 1 + α*Δt*θ accounts for implicit exchange at new time step
        !> This formulation ensures mass conservation and stability for θ ≥ 0.5
        conc_imm_new(j)=(conc_imm_old(j)*(1d0-this%imm_zones(j)%exch_rate*Delta_t*(1-theta))+this%imm_zones(j)%exch_rate*Delta_t*&
            (theta*conc_mob_new(j)+(1-theta)*conc_mob_old(j)))/(1+this%imm_zones(j)%exch_rate*Delta_t*theta) !< [M/L³] = ([M/L³]*[-] + [1/T]*[T]*[M/L³]) / [-]
    end do !< End loop over immobile zones
end subroutine !< End of compute_conc_imm_MRMT subroutine