!> @file params_aq_sol_m.f90
!> @brief Aqueous solution parameters module
!> @details This module defines the structure containing thermodynamic parameters 
!> for aqueous solutions, specifically Debye-Hückel parameters used in activity 
!> coefficient calculations
!> @author Generated documentation
!> @date 2025-11-18

!> @brief Aqueous solution parameters module
!> @details Contains the params_aq_sol_s type structure with Debye-Hückel parameters
module params_aq_sol_m
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module variables between invocations
    private !< Default visibility is private; expose only selected types/procedures
    !> @brief Aqueous solution parameters structure
    !> @details Structure containing Debye-Hückel equation parameters used for
    !> calculating activity coefficients in aqueous solutions. These parameters
    !> are temperature-dependent; default values are for 25°C (298.15 K)
    type, public :: params_aq_sol_s
        real(kind=8) :: A=0.5092 !< Debye-Hückel parameter A at 25°C (mol^-1/2·kg^1/2)
        real(kind=8) :: B=0.3283 !< Debye-Hückel parameter B at 25°C (mol^-1/2·kg^1/2·Å^-1)
    end type params_aq_sol_s !< End of aqueous solution parameters type definition
    
end module params_aq_sol_m !< End of aqueous solution parameters module