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
        real(kind=8), allocatable :: rk_new(:)              !< Kinetic reaction rates at end of time step
        real(kind=8), allocatable :: Rk_est(:)              !< Estimated kinetic reaction amounts
        real(kind=8), allocatable :: Rk(:)                  !< Kinetic reaction amount during a time step
        real(kind=8), allocatable :: Rk_accum(:)            !< Accumulated kinetic reaction amount during simulation
        real(kind=8), allocatable :: rk_mean(:)             !< Mean kinetic reaction rate during a time step
        real(kind=8), allocatable :: rk_old(:)              !< Kinetic reaction rates in previous time step
        real(kind=8), allocatable :: rk_old_old(:)          !< Kinetic reaction rates in two previous time steps
        real(kind=8), allocatable :: rk_old_old_old(:)      !< Kinetic reaction rates in three previous time steps
        real(kind=8), allocatable :: Re(:)                  !< Equilibrium reaction amount during a time step
        real(kind=8), allocatable :: re_mean(:)             !< Mean equilibrium reaction rate during a time step
        real(kind=8) :: volume=1d0                          !< Volume of solution or volumetric fraction (default 1) [L]
        integer(kind=4), allocatable :: var_act_species_indices(:) !< Indices of variable activity species
        integer(kind=4), allocatable :: cst_act_species_indices(:) !< Indices of constant activity species
    contains
        procedure :: set_name
        procedure :: set_id
        procedure :: set_density
        procedure :: set_pressure
        procedure :: set_temp
        procedure :: set_volume
        procedure :: set_conc_old
        procedure :: set_conc_old_old
        procedure :: set_re_mean
        procedure :: set_Rk_est
        procedure(set_concentrations), public, deferred :: set_concentrations
        procedure :: allocate_var_act_species_indices
        procedure :: allocate_cst_act_species_indices
        procedure :: update_rk_old
        procedure :: update_conc_old
        procedure :: copy_local_chemistry
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
            this%density = 0.9970749d0
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

    !> \brief Update kinetic reaction rates for previous time steps
    !> \param this Local chemistry object
    subroutine update_rk_old(this)
        class(local_chemistry_c) :: this
        this%rk_old_old_old = this%rk_old_old
        this%rk_old_old = this%rk_old
        this%rk_old = this%rk_new
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
        if (allocated(other%Rk)) then
            this%rk_new = other%rk_new
            this%Rk_est = other%Rk_est
            this%Rk = other%Rk
            this%Rk_accum = other%Rk_accum
            this%rk_mean = other%rk_mean
            this%rk_old = other%rk_old
            this%rk_old_old = other%rk_old_old
            this%rk_old_old_old = other%rk_old_old_old
        end if
        if (allocated(other%Re)) then
            this%Re = other%Re
            this%re_mean = other%re_mean
        end if
    end subroutine

    !> \brief Set the mean equilibrium reaction rates during a time step
    !> \param this Local chemistry object
    !> \param Delta_t Time step size
    subroutine set_re_mean(this,Delta_t)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in) :: Delta_t
        this%re_mean = this%Re / Delta_t
    end subroutine
    
    subroutine set_Rk_est(this,Rk_est)
        class(local_chemistry_c) :: this
        real(kind=8), intent(in) :: Rk_est(:)
        this%Rk_est = Rk_est
    end subroutine

end module