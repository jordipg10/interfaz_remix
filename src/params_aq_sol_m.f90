!> @file params_aq_sol_m.f90
!> @brief Aqueous solution parameters module
!> @details This module defines the structure containing thermodynamic parameters 
!> for aqueous solutions, specifically Debye-Hückel parameters used in activity 
!> coefficient calculations, and their evaluation at the temperature of a water
!> @author Generated documentation
!> @date 2025-11-18

!> @brief Aqueous solution parameters module
!> @details Contains the params_aq_sol_s type structure with Debye-Hückel parameters
module params_aq_sol_m
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module variables between invocations
    private !< Default visibility is private; expose only selected types/procedures
    public :: dielectric_cst_water
    !> @brief Aqueous solution parameters structure
    !> @details Structure containing Debye-Hückel equation parameters used for
    !> calculating activity coefficients in aqueous solutions. These parameters
    !> are temperature-dependent: the defaults are their values at 25°C (298.15 K) and
    !> 1 atm, and compute_params_aq_sol evaluates them at the temperature (and water
    !> density) of each water. They are therefore a property of each aqueous chemistry
    !> object (like its temperature), not of the shared aqueous phase.
    type, public :: params_aq_sol_s
        real(kind=8) :: A=0.510015d0 !< Debye-Hückel parameter A (default: compute_params_aq_sol at 25°C, 1 atm) [kg^1/2·mol^-1/2]
        real(kind=8) :: B=0.328489d0 !< Debye-Hückel parameter B (default: compute_params_aq_sol at 25°C, 1 atm) [kg^1/2·mol^-1/2·Å^-1]
    contains
        procedure :: compute_params_aq_sol !< Evaluate A and B at a temperature and water density
    end type params_aq_sol_s !< End of aqueous solution parameters type definition

    real(kind=8), parameter :: P_1atm_bar=1.01325d0 !< Pressure at which the dielectric constant is evaluated: 1 atm [bar]
    
contains
    !> @brief Relative dielectric constant (relative permittivity) of liquid water
    !> @details Bradley & Pitzer (1979), J. Phys. Chem. 83, 1599-1603 (the equation PHREEQC uses):
    !> @f[ \varepsilon(T,P) = D_{1000} + C \ln\frac{B+P}{B+1000}, \quad
    !>     D_{1000}=U_1 e^{U_2 T + U_3 T^2}, \quad C = U_4 + \frac{U_5}{U_6+T}, \quad B = U_7 + \frac{U_8}{T} + U_9 T @f]
    !> with T in K and P in bar. ε(25°C, 1 atm) = 78.38.
    !> @param[in] temp Temperature [K]
    !> @param[in] pressure Pressure [bar]
    !> @return eps Relative dielectric constant of water [-]
    pure function dielectric_cst_water(temp,pressure) result(eps)
        implicit none !< Enforce explicit variable declarations
        real(kind=8), intent(in) :: temp !< Temperature [K]
        real(kind=8), intent(in) :: pressure !< Pressure [bar]
        real(kind=8) :: eps !< Relative dielectric constant [-]
        real(kind=8), parameter :: U1=3.4279d2, U2=-5.0866d-3, U3=9.4690d-7, U4=-2.0525d0, U5=3.1159d3, &
            U6=-1.8289d2, U7=-8.0325d3, U8=4.2142d6, U9=2.1417d0 !< Bradley & Pitzer (1979) coefficients
        real(kind=8) :: D1000,C,Bp !< Terms of the Bradley-Pitzer equation
        D1000=U1*exp(temp*(U2+temp*U3)) !< Dielectric constant at 1000 bar
        C=U4+U5/(U6+temp)
        Bp=U7+U8/temp+U9*temp
        eps=D1000+C*log((Bp+pressure)/(Bp+1d3))
    end function

    !> @brief Evaluate the Debye-Hückel parameters A and B at a temperature
    !> @details PHREEQC convention (Parkhurst & Appelo, 2013):
    !> @f[ A = 1.82483\cdot10^{6}\,\frac{\sqrt{\rho_w}}{(\varepsilon T)^{3/2}}, \qquad
    !>     B = 50.2916\,\frac{\sqrt{\rho_w}}{(\varepsilon T)^{1/2}} @f]
    !> with T the temperature [K], ρ_w the density of pure water [kg/L] and ε its relative
    !> dielectric constant (Bradley & Pitzer, 1979, at 1 atm). At 25°C (ρ_w = 0.99705 kg/L)
    !> this gives A = 0.5100 and B = 0.3285 (the defaults of the type); at 60°C A = 0.5459
    !> and B = 0.3345.
    !> @param[inout] this Aqueous solution parameters
    !> @param[in] temp Temperature [K]
    !> @param[in] density Density of pure water at that temperature [kg/L]
    subroutine compute_params_aq_sol(this,temp,density)
        implicit none !< Enforce explicit variable declarations
        class(params_aq_sol_s), intent(inout) :: this !< Aqueous solution parameters [-]
        real(kind=8), intent(in) :: temp !< Temperature [K]
        real(kind=8), intent(in) :: density !< Density of pure water [kg/L]
        real(kind=8) :: eps_T !< Product of the dielectric constant and the temperature [K]
        eps_T=dielectric_cst_water(temp,P_1atm_bar)*temp
        this%A=1.82483d6*sqrt(density)/(eps_T*sqrt(eps_T)) !< [kg^1/2·mol^-1/2]
        this%B=50.2916d0*sqrt(density)/sqrt(eps_T) !< [kg^1/2·mol^-1/2·Å^-1]
    end subroutine
    
end module params_aq_sol_m !< End of aqueous solution parameters module