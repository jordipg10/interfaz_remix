!> Local chemistry abstract superclass
!!> This class contains the local chemical information
module local_chemistry_m
    implicit none
    save
    type, public, abstract :: local_chemistry_c
        character(len=256) :: name !> name of local chemistry
        integer(kind=4) :: id !> id of local chemistry (for labeling purposes) (must be non-negative)
        real(kind=8) :: temp !> temperature [K]
        real(kind=8) :: density !> [kg/L]
        real(kind=8) :: pressure !> [atm]
        real(kind=8), allocatable :: concentrations(:) !> species concentrations
        real(kind=8), allocatable :: conc_old(:) !> species concentrations at previous time step
        real(kind=8), allocatable :: conc_old_old(:) !> species concentrations at two previous time steps (chapuza)
        real(kind=8), allocatable :: activities(:) !> species activities
        real(kind=8), allocatable :: log_act_coeffs(:) !> log_10 activity coefficients
        real(kind=8), allocatable :: log_Jacobian_act_coeffs(:,:) !> log_10 Jacobian activity coefficients
        real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
        real(kind=8), allocatable :: Rk_est(:) !> estimated kinetic reaction amounts
        real(kind=8), allocatable :: Rk_mean(:) !> mean kinetic reaction amount during a time step
        real(kind=8), allocatable :: Re_mean(:) !> mean equilibrium reaction amount during a time step
        real(kind=8), allocatable :: rk_old(:) !> kinetic reaction rates in previous time step
        real(kind=8), allocatable :: rk_old_old(:) !> kinetic reaction rates in two previous time step (chapuza)
        real(kind=8), allocatable :: rk_old_old_old(:) !> kinetic reaction rates in three previous time step (chapuza)
        real(kind=8), allocatable :: r_eq(:) !> equilibrium reaction rates
        real(kind=8) :: volume=1d0 !> volume of solution or volumetric fraction (1 by default)
        integer(kind=4), allocatable :: var_act_species_indices(:)
        integer(kind=4), allocatable :: cst_act_species_indices(:)
    contains
    !> Set
        procedure, public :: set_name
        procedure, public :: set_id
        procedure, public :: set_density
        procedure, public :: set_pressure
        procedure, public :: set_temp
        procedure, public :: set_volume
        procedure, public :: set_conc_old
        procedure, public :: set_conc_old_old
        procedure(set_concentrations), public, deferred :: set_concentrations
    !> Allocate
        procedure, public :: allocate_var_act_species_indices
        procedure, public :: allocate_cst_act_species_indices
        !procedure, public :: allocate_reaction_rates
    !> Update
        procedure, public :: update_rk_old
        procedure, public :: update_conc_old
    end type
    
    abstract interface
        subroutine set_concentrations(this,conc)
            import local_chemistry_c
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
        end subroutine
    end interface
        
    contains
        subroutine set_name(this,name)
        !> Set name of local chemistry
            implicit none
            class(local_chemistry_c) :: this
            character(len=*), intent(in) :: name !> name of local chemistry
            this%name=name
        end subroutine

        subroutine set_temp(this,temp)
        !> Set temperature of local chemistry
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: temp !> Kelvin
            if (present(temp)) then
                this%temp=temp
            else
                this%temp=298.15 !> default
            end if
        end subroutine
        
        subroutine set_volume(this,vol)
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: vol !> [L]
            if (present(vol)) then
                this%volume=vol
            else
                this%volume=1d0 !> default
            end if
        end subroutine
        
        subroutine set_pressure(this,pressure)
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: pressure
            if (present(pressure)) then
                this%pressure=pressure
            else
                this%pressure=1d0 !> default
            end if
        end subroutine
        
        subroutine set_density(this,density)
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: density
            if (present(density)) then
                this%density=density
            else
                this%density=0.9970749 !> [kg/L] (water)
            end if
        end subroutine
        
        !subroutine set_concentrations(this,conc)
        !    implicit none
        !    class(local_chemistry_c) :: this
        !    real(kind=8), intent(in) :: conc(:)
        !    this%concentrations=conc
        !end subroutine
        
        subroutine allocate_var_act_species_indices(this,num_var)
            implicit none
            class(local_chemistry_c) :: this
            integer(kind=4), intent(in) :: num_var
            allocate(this%var_act_species_indices(num_var))
        end subroutine
        
        subroutine allocate_cst_act_species_indices(this,num_cst)
            implicit none
            class(local_chemistry_c) :: this
            integer(kind=4), intent(in) :: num_cst
            allocate(this%cst_act_species_indices(num_cst))
        end subroutine
        
        subroutine update_rk_old(this)
        implicit none
        class(local_chemistry_c) :: this
        this%rk_old_old_old=this%rk_old_old
        this%rk_old_old=this%rk_old
        this%rk_old=this%rk
        end subroutine

        subroutine update_conc_old(this)
            implicit none
            class(local_chemistry_c) :: this
            if (allocated(this%concentrations)) then
                this%conc_old_old=this%conc_old
                this%conc_old=this%concentrations
            end if
        end subroutine
        
        subroutine set_conc_old(this,conc_old)
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: conc_old(:)
            !this%conc_old_old=this%conc_old
            if (present(conc_old)) then
                this%conc_old=conc_old
            else if (allocated(this%concentrations)) then
                this%conc_old=this%concentrations
            end if
        end subroutine
        
        subroutine set_conc_old_old(this,conc_old_old)
            implicit none
            class(local_chemistry_c) :: this
            real(kind=8), intent(in), optional :: conc_old_old(:)
            !this%conc_old_old=this%conc_old
            if (present(conc_old_old)) then
                this%conc_old_old=conc_old_old
            else if (allocated(this%conc_old)) then
                this%conc_old_old=this%conc_old
            end if
        end subroutine

        ! subroutine allocate_reaction_rates(this)
        !     implicit none
        !     class(local_chemistry_c) :: this
        !     !integer(kind=4), intent(in) :: num_reactions
        !     ! allocate(this%rk(num_reactions))
        !     ! allocate(this%Rk_est(num_reactions))
        !     ! allocate(this%Rk_mean(num_reactions))
        !     ! allocate(this%Re_mean(num_reactions))
        !     ! allocate(this%rk_old(num_reactions))
        !     ! allocate(this%rk_old_old(num_reactions))
        !     ! allocate(this%rk_old_old_old(num_reactions))
        !     ! allocate(this%r_eq(num_reactions))
        ! end subroutine
        
        subroutine set_id(this,id)
        implicit none
        class(local_chemistry_c) :: this
        integer(kind=4), intent(in) :: id !> id of local chemistry
        if (id<0) error stop "ID must be non-negative"
        this%id=id
        end subroutine
end module