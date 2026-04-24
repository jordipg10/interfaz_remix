!> @file aq_species_m.f90
!> @brief Aqueous species subclass module
!> @details This module defines the aqueous species class which extends the base species class.
!>          It contains properties specific to aqueous species and implements various models
!>          for computing logarithmic activity coefficients.
!> @author [Author name]
!> @date [Date]

!> @brief Aqueous species module
!> @details Contains the properties of aqueous species and computes logarithm activity 
!>          coefficient with different models (Debye-Hückel, Davies, Truesdell-Jones)
module aq_species_m
    use species_m, only : species_c                                   !< Import parent species module for inheritance
    implicit none                                   !< Require explicit variable declarations
    save                                           !< Preserve module variables between procedure calls
    private                                        !< Default accessibility is private
    !> @brief Aqueous species derived type
    !> @details Extends the base species_c class with properties and methods specific to
    !>          aqueous species, including activity coefficient calculations
    type, public, extends(species_c) :: aq_species_c
        real(kind=8) :: ionic_radius               !< Ionic radius expressed in angstroms [Å]
        real(kind=8) :: alk_contrib                !< Contribution to total alkalinity
    contains
        !> @name Activity coefficient calculation methods
        !> @{
        ! procedure :: Debye_Huckel_restr    !< Restricted Debye-Hückel model
        ! procedure :: Debye_Huckel_ampl     !< Extended Debye-Hückel model  
        ! procedure :: Davies                !< Davies equation model
        ! procedure :: Truesdell_Jones      !< Truesdell-Jones model
        !> @}
    end type
    contains !< Module procedures implementation section
        
        !> @brief Restricted Debye-Hückel activity coefficient model
        !> @details Implements the restricted Debye-Hückel equation for calculating
        !>          activity coefficients: log γ = -A·√I·z²
        !>          Valid for ionic strengths < 0.1 M
        !> @param[in] this Aqueous species object
        !> @param[in] ionic_strength Ionic strength of the solution [M]
        !> @param[in] A Debye-Hückel constant A [M^(-1/2)]
        !> @param[in] B Debye-Hückel constant B [M^(-1/2) Å^(-1)]
        !> @param[out] log_act_coeff Base-10 logarithm of activity coefficient [-]
        ! subroutine Debye_Huckel_restr(this,ionic_strength,A,B,log_act_coeff)
        !     implicit none
        !     class(aq_species_c) :: this                              !< Aqueous species object
        !     real(kind=8), intent(in) :: ionic_strength               !< Ionic strength [M]
        !     real(kind=8), intent(in) :: A,B                          !< Debye-Hückel constants
        !     real(kind=8), intent(out) :: log_act_coeff               !< log₁₀(γ) [-]
            
        !     !> Calculate restricted Debye-Hückel activity coefficient
        !     log_act_coeff=-A*sqrt(ionic_strength)*this%valence**2
        ! end subroutine
        
        !> @brief Extended Debye-Hückel activity coefficient model
        !> @details Implements the extended Debye-Hückel equation with ion size parameter:
        !>          log γ = -A·√I·z² / (1 + a·B·√I)
        !>          Valid for ionic strengths up to ~0.5 M
        !> @param[in] this Aqueous species object
        !> @param[in] ionic_strength Ionic strength of the solution [M]
        !> @param[in] A Debye-Hückel constant A [M^(-1/2)]
        !> @param[in] B Debye-Hückel constant B [M^(-1/2) Å^(-1)]
        !> @param[out] log_act_coeff Base-10 logarithm of activity coefficient [-]
        ! subroutine Debye_Huckel_ampl(this,ionic_strength,A,B,log_act_coeff)
        !     implicit none
        !     class(aq_species_c) :: this                              !< Aqueous species object
        !     real(kind=8), intent(in) :: ionic_strength               !< Ionic strength [M]
        !     real(kind=8), intent(in) :: A,B                          !< Debye-Hückel constants
        !     real(kind=8), intent(out) :: log_act_coeff               !< log₁₀(γ) [-]
            
        !     !> Calculate extended Debye-Hückel activity coefficient with ion size correction
        !     log_act_coeff=-(A*sqrt(ionic_strength)*this%valence**2)/(&
        !         1d0+this%params_act_coeff%ion_size_param*B*sqrt(ionic_strength))
        ! end subroutine
        
        !> @brief Davies equation activity coefficient model
        !> @details Implements the Davies equation, a semi-empirical modification of
        !>          Debye-Hückel: log γ = -A·z²·[(√I/(1+√I)) - 0.3·I]
        !>          Valid for intermediate ionic strengths (0.1-0.5 M)
        !> @param[in] this Aqueous species object
        !> @param[in] ionic_strength Ionic strength of the solution [M]
        !> @param[in] A Debye-Hückel constant A [M^(-1/2)]
        !> @param[in] B Debye-Hückel constant B (not used in Davies equation)
        !> @param[out] log_act_coeff Base-10 logarithm of activity coefficient [-]
        ! subroutine Davies(this,ionic_strength,A,B,log_act_coeff)
        !     implicit none
        !     class(aq_species_c) :: this                              !< Aqueous species object
        !     real(kind=8), intent(in) :: ionic_strength               !< Ionic strength [M]
        !     real(kind=8), intent(in) :: A,B                          !< Debye-Hückel constants (B not used)
        !     real(kind=8), intent(out) :: log_act_coeff               !< log₁₀(γ) [-]
            
        !     !> Calculate Davies activity coefficient
        !     log_act_coeff=-A*(this%valence**2)*((sqrt(ionic_strength)/(1d0+sqrt(ionic_strength)))-3d-1*ionic_strength)
        ! end subroutine
        
        !> @brief Truesdell-Jones activity coefficient model
        !> @details Implements the Truesdell-Jones equation with species-specific parameters:
        !>          log γ = -A·√I·z²/(1+B·a_TJ·√I) + b_TJ·I
        !>          Uses species-specific a_TJ and b_TJ coefficients for improved accuracy
        !> @param[in] this Aqueous species object containing a_TJ and b_TJ parameters
        !> @param[in] ionic_strength Ionic strength of the solution [M]
        !> @param[in] A Debye-Hückel constant A [M^(-1/2)]
        !> @param[in] B Debye-Hückel constant B [M^(-1/2) Å^(-1)]
        !> @param[out] log_act_coeff Base-10 logarithm of activity coefficient [-]
        ! subroutine Truesdell_Jones(this,ionic_strength,A,B,log_act_coeff)
        !     implicit none
        !     class(aq_species_c) :: this                              !< Aqueous species object
        !     real(kind=8), intent(in) :: ionic_strength               !< Ionic strength [M]
        !     real(kind=8), intent(in) :: A,B                          !< Debye-Hückel constants
        !     real(kind=8), intent(out) :: log_act_coeff               !< log₁₀(γ) [-]
            
        !     !> Calculate Truesdell-Jones activity coefficient using species-specific parameters
        !     log_act_coeff=-A*sqrt(ionic_strength)*(this%valence**2)/(1d0+B*this%params_act_coeff%a_TJ*sqrt(ionic_strength))+&
        !     this%params_act_coeff%b_TJ*ionic_strength
        ! end subroutine
    
    
end module !< End of aq_species_m module