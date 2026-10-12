!> \file local_chemistry_m.f90
!> \brief Defines the abstract superclass for local chemical information in simulations.
!> \details
!>   - Provides the abstract type \ref local_chemistry_c for storing local chemical state and properties.
!>   - Includes procedures for setting, allocating, updating, and assigning chemical data.
!>   - Used as a base for concrete local chemistry implementations.
!> Author: Jordi Petchamé-Guerrero
!> Date: 2025
module local_chemistry_m
    implicit none
    save
    private
    !> \brief Abstract superclass for local chemical information.
    type, public, abstract :: local_chemistry_c
        character(len=256) :: name         !< Name of local chemistry
        integer(kind=4) :: id              !< ID of local chemistry (non-negative)
        real(kind=8) :: temp               !< Temperature [K]
        real(kind=8) :: density            !< Density
        real(kind=8) :: pressure           !< Pressure
        real(kind=8), allocatable :: concentrations(:)      !< Species concentrations at current time step
        real(kind=8), allocatable :: conc_old(:)            !< Species concentrations after previous time step
        real(kind=8), allocatable :: conc_old_old(:)        !< Species concentrations before previous time step
        real(kind=8), allocatable :: activities(:)          !< Species activities
        real(kind=8), allocatable :: log_act_coeffs(:)      !< log_10 activity coefficients
        real(kind=8), allocatable :: log_Jacobian_act_coeffs(:,:) !< log_10 Jacobian activity coefficients
        real(kind=8), allocatable :: rxn_rate_kin(:)        !< Kinetic reaction rates at end of time step
        real(kind=8), allocatable :: rxn_amnt_kin(:)        !< Kinetic reaction amount during a time step
        real(kind=8), allocatable :: cum_rxn_amnt_kin(:)      !< Cumulative kinetic reaction amount over the simulation
        real(kind=8), allocatable :: cum_rxn_amnt_eq(:)       !< Cumulative equilibrium reaction amount over the simulation
        real(kind=8), allocatable :: rxn_rate_kin_avg(:)      !< Mean kinetic reaction rate during a time step
        real(kind=8), allocatable :: rxn_rate_kin_old(:)      !< Kinetic reaction rates after previous time step
        real(kind=8), allocatable :: rk_old_old(:)          !< Kinetic reaction rates in two previous time steps
        real(kind=8), allocatable :: rk_old_old_old(:)      !< Kinetic reaction rates in three previous time steps
        real(kind=8), allocatable :: rxn_amnt_eq(:)         !< Equilibrium reaction amount during a time step
        real(kind=8), allocatable :: rxn_rate_eq(:)         !< equilibrium reaction rates during a time step
        real(kind=8), allocatable :: rxn_rate_eq_old(:)      !< equilibrium reaction rates in previous time step
        real(kind=8), allocatable :: rxn_rate_eq_old_old(:)      !< equilibrium reaction rates in two previous time steps
        real(kind=8), allocatable :: rxn_rate_eq_old_old_old(:)  !< equilibrium reaction rates in three previous time steps
        real(kind=8), allocatable :: contrib_rxn_amnt_eq(:)      !< Contribution of equilibrium reaction amounts to species mass balance (rho_eq); NOT the reaction amount Re itself
        real(kind=8), allocatable :: contrib_rxn_amnt_kin(:)     !< Contribution of kinetic reaction amounts to species mass balance (rho_kin); NOT the reaction amount Rk itself
        real(kind=8) :: volume=1d0                          !< Volume of solution or volumetric fraction (default 1) [L]
        integer(kind=4), allocatable :: var_act_species_indices(:) !< Indices of variable activity species in chemical system
        integer(kind=4), allocatable :: cst_act_species_indices(:) !< Indices of constant activity species in chemical system
    contains
        procedure :: set_name
        procedure :: set_id
        procedure :: set_density
        procedure :: set_pressure
        procedure :: set_temp
        procedure :: set_volume
        procedure :: set_conc_old
        procedure :: set_conc_old_old
        ! procedure :: set_eq_rxn_rates
        procedure(set_concentrations), public, deferred :: set_concentrations
        procedure :: allocate_var_act_species_indices
        procedure :: allocate_cst_act_species_indices
        procedure :: update_rxn_rates_and_amounts
        procedure :: update_cum_rxn_amnt
        procedure :: update_conc_old
        procedure :: copy_local_chemistry
        ! procedure :: compute_rxn_rate_kin_avg
    end type local_chemistry_c
    
    abstract interface
        !> \brief Abstract deferred procedure to set concentrations
        !> \param this Local chemistry object
        !> \param conc Concentrations array
        subroutine set_concentrations(this,conc)
            import local_chemistry_c
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
        end subroutine
    end interface
        
contains

    !> \brief Set the name of local chemistry
    !> \param this Local chemistry object
    !> \param name Name to assign
    subroutine set_name(this,name)
        class(local_chemistry_c) :: this
        character(len=*), intent(in) :: name
        this%name = name
    end subroutine

    !> \brief Set the temperature of local chemistry
    !> \param this Local chemistry object
    !> \param temp [optional] Temperature in Kelvin
    subroutine set_temp(this,temp)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: temp
        if (present(temp)) then
            this%temp = temp
        else
            this%temp = 298.15d0
        end if
    end subroutine

    !> \brief Set the volume of local chemistry
    !> \param this Local chemistry object
    !> \param vol [optional] Volume in m^3
    subroutine set_volume(this,vol)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: vol
        if (present(vol)) then
            this%volume = vol
        else
            this%volume = 1d0
        end if
    end subroutine

    !> \brief Set the pressure of local chemistry
    !> \param this Local chemistry object
    !> \param pressure [optional] Pressure in atm
    subroutine set_pressure(this,pressure)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: pressure
        if (present(pressure)) then
            this%pressure = pressure
        else
            this%pressure = 1d0
        end if
    end subroutine

    !> \brief Set the density of local chemistry
    !> \param this Local chemistry object
    !> \param density [optional] Density in kg/L
    subroutine set_density(this,density)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: density
        if (present(density)) then
            this%density = density
        else
            this%density = 0.9970749d0 !< Default density of water at 25°C [kg/L]
        end if
    end subroutine

    !> \brief Allocate indices for variable activity species
    !> \param this Local chemistry object
    !> \param num_var Number of variable activity species
    subroutine allocate_var_act_species_indices(this,num_var)
        class(local_chemistry_c) :: this
        integer(kind=4), intent(in) :: num_var
        allocate(this%var_act_species_indices(num_var))
    end subroutine

    !> \brief Allocate indices for constant activity species
    !> \param this Local chemistry object
    !> \param num_cst Number of constant activity species
    subroutine allocate_cst_act_species_indices(this,num_cst)
        class(local_chemistry_c) :: this
        integer(kind=4), intent(in) :: num_cst
        allocate(this%cst_act_species_indices(num_cst))
    end subroutine

    !> \brief Update kinetic and equilibrium reaction-rate histories (shift current rates into the old-step arrays)
    !> \param this Local chemistry object
    subroutine update_rxn_rates_and_amounts(this)
        class(local_chemistry_c) :: this
        !> Guard the kinetic-rate history shift: dummy solid chemistries (e.g. gas-only waters
        !> with no minerals/kinetics) leave these arrays unallocated.
        if (allocated(this%rk_old_old)) this%rk_old_old_old = this%rk_old_old
        if (allocated(this%rxn_rate_kin_old)) this%rk_old_old = this%rxn_rate_kin_old
        if (allocated(this%rxn_rate_kin)) this%rxn_rate_kin_old = this%rxn_rate_kin
        !if (allocated(this%rxn_amnt_kin)) this%rxn_rate_kin_avg = 0.5d0*(this%rxn_rate_kin + this%rxn_rate_kin_old)
        !if (allocated(this%rxn_amnt_kin)) this%cum_rxn_amnt_kin = this%cum_rxn_amnt_kin + this%rxn_amnt_kin
        !> Shift the equilibrium mean-rate history end-of-step, in lockstep with the kinetic history,
        !> so that get_r_old/get_r_old_old/get_r_old_old_old pair rxn_rate_eq_old/rxn_rate_eq_old_old/rxn_rate_eq_old_old_old
        !> with rk_old/rk_old_old/rk_old_old_old at the same time level (re_mean itself is the current
        !> value, set per water in set_eq_rxn_rates). Guard each shift by its source: a dummy/gas-only solid
        !> may have re_mean allocated without the _old history arrays.
        if (allocated(this%rxn_rate_eq_old_old)) this%rxn_rate_eq_old_old_old = this%rxn_rate_eq_old_old
        if (allocated(this%rxn_rate_eq_old)) this%rxn_rate_eq_old_old = this%rxn_rate_eq_old
        if (allocated(this%rxn_rate_eq)) this%rxn_rate_eq_old = this%rxn_rate_eq
        !if (allocated(this%rxn_amnt_eq)) this%cum_rxn_amnt_eq = this%cum_rxn_amnt_eq + this%rxn_amnt_eq
    end subroutine

    !> \brief Update cumulative reaction amounts for kinetic and equilibrium reactions
    !> \param this Local chemistry object
    subroutine update_cum_rxn_amnt(this) 
        class(local_chemistry_c) :: this
        !> Guard the kinetic-rate history shift: dummy solid chemistries (e.g. gas-only waters
        !> with no minerals/kinetics) leave these arrays unallocated.
        ! if (allocated(this%rk_old_old)) this%rk_old_old_old = this%rk_old_old
        ! if (allocated(this%rxn_rate_kin_old)) this%rk_old_old = this%rxn_rate_kin_old
        ! if (allocated(this%rxn_rate_kin)) this%rxn_rate_kin_old = this%rxn_rate_kin
        !if (allocated(this%rxn_amnt_kin)) this%rxn_rate_kin_avg = 0.5d0*(this%rxn_rate_kin + this%rxn_rate_kin_old)
        if (allocated(this%cum_rxn_amnt_kin)) this%cum_rxn_amnt_kin = this%cum_rxn_amnt_kin + this%rxn_amnt_kin
        !> Shift the equilibrium mean-rate history end-of-step, in lockstep with the kinetic history,
        !> so that get_r_old/get_r_old_old/get_r_old_old_old pair rxn_rate_eq_old/rxn_rate_eq_old_old/rxn_rate_eq_old_old_old
        !> with rk_old/rk_old_old/rk_old_old_old at the same time level (re_mean itself is the current
        !> value, set per water in set_eq_rxn_rates). Guard each shift by its source: a dummy/gas-only solid
        !> may have re_mean allocated without the _old history arrays.
        ! if (allocated(this%rxn_rate_eq_old_old)) this%rxn_rate_eq_old_old_old = this%rxn_rate_eq_old_old
        ! if (allocated(this%rxn_rate_eq_old)) this%rxn_rate_eq_old_old = this%rxn_rate_eq_old
        ! if (allocated(this%rxn_rate_eq)) this%rxn_rate_eq_old = this%rxn_rate_eq
        if (allocated(this%cum_rxn_amnt_eq)) this%cum_rxn_amnt_eq = this%cum_rxn_amnt_eq + this%rxn_amnt_eq
    end subroutine

    !> \brief Update concentrations for previous time steps
    !> \param this Local chemistry object
    subroutine update_conc_old(this)
        class(local_chemistry_c) :: this
        if (allocated(this%concentrations)) then
            this%conc_old_old = this%conc_old
            this%conc_old = this%concentrations
        end if
    end subroutine

    !> \brief Set concentrations at previous time step
    !> \param this Local chemistry object
    !> \param conc_old [optional] Concentrations array
    subroutine set_conc_old(this,conc_old)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: conc_old(:)
        if (present(conc_old)) then
            this%conc_old = conc_old
        else if (allocated(this%concentrations)) then
            this%conc_old = this%concentrations
        end if
    end subroutine

    !> \brief Set concentrations at two previous time steps
    !> \param this Local chemistry object
    !> \param conc_old_old [optional] Concentrations array
    subroutine set_conc_old_old(this,conc_old_old)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in), optional :: conc_old_old(:)
        if (present(conc_old_old)) then
            this%conc_old_old = conc_old_old
        else if (allocated(this%conc_old)) then
            this%conc_old_old = this%conc_old
        end if
    end subroutine

    !> \brief Set the ID of local chemistry
    !> \param this Local chemistry object
    !> \param id ID to assign (must be non-negative)
    !> \throws error stop if id < 0
    subroutine set_id(this,id)
        class(local_chemistry_c) :: this
        integer(kind=4), intent(in) :: id
        if (id < 0) error stop "ID must be non-negative"
        this%id = id
    end subroutine

    !> \brief Assign all local chemistry properties from another object
    !> \param this Local chemistry object
    !> \param other Source local chemistry object
    subroutine copy_local_chemistry(this,other)
        class(local_chemistry_c) :: this
        class(local_chemistry_c), intent(in) :: other
        this%id = other%id
        this%name = other%name
        this%pressure = other%pressure
        this%density = other%density
        this%temp = other%temp
        this%volume = other%volume
        if (allocated(other%concentrations)) then
            this%concentrations = other%concentrations
            this%conc_old = other%conc_old
            this%conc_old_old = other%conc_old_old
            !this%contrib_rxn_amnt_kin = other%contrib_rxn_amnt_kin
            !this%contrib_rxn_amnt_eq = other%contrib_rxn_amnt_eq
        end if
        if (allocated(other%contrib_rxn_amnt_kin)) then
            this%contrib_rxn_amnt_kin = other%contrib_rxn_amnt_kin
        end if
        if (allocated(other%contrib_rxn_amnt_eq)) then
            this%contrib_rxn_amnt_eq = other%contrib_rxn_amnt_eq
        end if
        if (allocated(other%activities)) then
            this%activities = other%activities
        end if
        if (allocated(other%log_act_coeffs) .and. allocated(other%log_Jacobian_act_coeffs)) then
            this%log_act_coeffs = other%log_act_coeffs
            this%log_Jacobian_act_coeffs = other%log_Jacobian_act_coeffs
        end if
        if (allocated(other%var_act_species_indices) .and. allocated(other%cst_act_species_indices)) then
            this%var_act_species_indices = other%var_act_species_indices
            this%cst_act_species_indices = other%cst_act_species_indices
        end if
        if (allocated(other%rxn_amnt_kin)) then
            this%rxn_rate_kin = other%rxn_rate_kin
            this%rxn_amnt_kin = other%rxn_amnt_kin
            this%cum_rxn_amnt_kin = other%cum_rxn_amnt_kin
            this%rxn_rate_kin_avg = other%rxn_rate_kin_avg
            this%rxn_rate_kin_old = other%rxn_rate_kin_old
            this%rk_old_old = other%rk_old_old
            this%rk_old_old_old = other%rk_old_old_old
            !this%contrib_rxn_amnt_kin = other%contrib_rxn_amnt_kin
        end if
        !> the eq arrays are not always allocated together (solid/gas allocators vs the aqueous one),
        !> so copy each on its own: assigning from an unallocated allocatable is undefined behaviour
        if (allocated(other%rxn_amnt_eq)) then
            this%rxn_amnt_eq = other%rxn_amnt_eq
        end if
        if (allocated(other%rxn_rate_eq)) then
            this%rxn_rate_eq = other%rxn_rate_eq
        end if
        if (allocated(other%cum_rxn_amnt_eq)) then
            this%cum_rxn_amnt_eq = other%cum_rxn_amnt_eq
        end if
        if (allocated(other%rxn_rate_eq_old)) then
            this%rxn_rate_eq_old = other%rxn_rate_eq_old
        end if
        if (allocated(other%rxn_rate_eq_old_old)) then
            this%rxn_rate_eq_old_old = other%rxn_rate_eq_old_old
        end if
        if (allocated(other%rxn_rate_eq_old_old_old)) then
            this%rxn_rate_eq_old_old_old = other%rxn_rate_eq_old_old_old
        end if
    end subroutine

end module