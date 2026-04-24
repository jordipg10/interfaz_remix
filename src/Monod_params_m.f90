!> \file Monod_params_m.f90
!> \brief Monod kinetic parameters module for microbial reaction rates
!> \details
!>   Defines parameters for Monod kinetic rate laws used in microbial reactions.
!>   
!>   Monod kinetics describe substrate-limited microbial growth:
!>   - Rate = μ_max * [S]/(K_M + [S]) * [Biomass]
!>   - Includes half-saturation constants (K_M)
!>   - Handles multiple electron donors/acceptors
!>   - Supports inhibition terms
!>   
!>   Applications:
!>   - Aerobic/anaerobic respiration
!>   - Fermentation reactions
!>   - Biodegradation processes
!>   - Denitrification, sulfate reduction, etc.
!>   
!>   Extends kin_params_c to inherit general kinetic parameter functionality.
!>
!> \author jordi Petchamé-Guerrero
!> \date October 2025

module Monod_params_m
    use kin_params_m, only: kin_params_c                                      !< Base kinetic parameters class
    implicit none
    save
    
    !> \brief Monod kinetic parameters class
    !> \details
    !>   Contains parameters for Monod-type microbial kinetic rate expressions.
    !>   
    !>   Monod rate law structure:
    !>   \f[
    !>   r = k \cdot \frac{[S_1]}{K_{M,1} + [S_1]} \cdot \frac{[S_2]}{K_{M,2} + [S_2]} 
    !>       \cdot \prod_{i} \frac{K_{inh,i}}{K_{inh,i} + [I_i]} \cdot [Biomass]
    !>   \f]
    !>   
    !>   Where:
    !>   - k: maximum specific growth rate (inherited from kin_params_c)
    !>   - S₁, S₂: electron donor and acceptor concentrations
    !>   - K_M: half-saturation constants (Monod constants)
    !>   - I_i: inhibitor concentrations
    !>   - K_inh: inhibition constants
    !>   
    !>   Number of terms:
    !>   - Always 2 substrate terms (donor + acceptor)
    !>   - Variable number of inhibition terms (n_inh)
    !>   - Total terms = n_inh + 2
    !>   
    !>   Threshold concentrations:
    !>   - Minimum concentrations below which reaction doesn't occur
    !>   - Prevents numerical issues at very low concentrations
    type, public, extends(kin_params_c) :: Monod_params_c
        integer(kind=4) :: num_terms                                          !< Total number of terms (2 substrates + n_inh inhibitors)
        integer(kind=4) :: n_inh                                              !< Number of inhibitor species
        !real(kind=8), allocatable :: conc_thr(:)                              !< Threshold concentrations for electron donors & acceptors
        real(kind=8), allocatable :: k_M(:)                                   !< Half-saturation (Monod) constants for substrates [mol/L]
        real(kind=8), allocatable :: k_inh(:)                                 !< Inhibition constants [mol/L]
        !real(kind=8), allocatable :: conc_thr_inh(:)                          !< Threshold concentrations for inhibitors
        
    contains
        procedure :: allocate_k_inh                                   !< Allocate inhibition constants array
        procedure :: allocate_k_M                                     !< Allocate Monod constants array
        procedure :: compute_num_terms                                !< Compute total number of terms
    end type
        
    !> \brief Interface declarations (currently empty)
    !> \details
    !>   Reserved for future external interface definitions.
    !>   Currently no external interfaces needed for Monod parameters.
    interface
     
    end interface
    
    contains
      
        !> \brief Allocate Monod (half-saturation) constants array
        !> \details
        !>   Allocates k_M array for two substrate terms:
        !>   - k_M(1): half-saturation constant for electron donor
        !>   - k_M(2): half-saturation constant for electron acceptor
        !>   
        !>   The Monod term: [S]/(K_M + [S])
        !>   - Approaches 0 when [S] << K_M (substrate-limited)
        !>   - Approaches 1 when [S] >> K_M (substrate-saturated)
        !>   - Equals 0.5 when [S] = K_M (half-maximum rate)
        !>
        !> \param[inout] this Monod parameters object
        subroutine allocate_k_M(this)
            implicit none
            class(Monod_params_c) :: this                                     !< Monod parameters object
            allocate(this%k_M(2))                                             !< Always 2 substrates (donor + acceptor)
        end subroutine
        
        !> \brief Allocate inhibition constants array
        !> \details
        !>   Allocates k_inh array for inhibitor species.
        !>   Optionally sets number of inhibitors before allocation.
        !>   
        !>   Inhibition term: K_inh/(K_inh + [I])
        !>   - Approaches 1 when [I] << K_inh (no inhibition)
        !>   - Approaches 0 when [I] >> K_inh (strong inhibition)
        !>   - Equals 0.5 when [I] = K_inh (half-maximum inhibition)
        !>   
        !>   Common inhibitors:
        !>   - Oxygen for anaerobic processes
        !>   - Product inhibition (e.g., ammonia in nitrification)
        !>   - Toxic substances
        !>
        !> \param[inout] this Monod parameters object
        !> \param[in] n_inh Number of inhibitors (optional, must be >= 0)
        subroutine allocate_k_inh(this,n_inh)
            implicit none
            class(Monod_params_c) :: this                                     !< Monod parameters object
            integer(kind=4), intent(in), optional :: n_inh                    !< Number of inhibitors
            if (present(n_inh)) then
                if (n_inh<0) then                                             !< Validate non-negative
                    error stop "Number of inhibitors must be positive"
                else
                    this%n_inh=n_inh                                          !< Set inhibitor count
                end if
            end if
            allocate(this%k_inh(this%n_inh))                                  !< Allocate inhibition constants array
        end subroutine
        
        !> \brief Compute total number of Monod terms
        !> \details
        !>   Calculates total terms in Monod rate expression:
        !>   - Always 2 substrate terms (electron donor and acceptor)
        !>   - Plus n_inh inhibition terms
        !>   
        !>   Total terms = 2 (substrates) + n_inh (inhibitors)
        !>   
        !>   Used to:
        !>   - Size arrays for term evaluation
        !>   - Allocate concentration arrays
        !>   - Validate input data consistency
        !>
        !> \param[inout] this Monod parameters object
        subroutine compute_num_terms(this)
            implicit none
            class(Monod_params_c) :: this                                     !< Monod parameters object
            this%num_terms=this%n_inh+2                                       !< Total = inhibitors + 2 substrates
        end subroutine
        
end module