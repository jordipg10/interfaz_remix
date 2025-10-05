!> Gas subclass: contains attributes gases
module gas_m
    use species_m
    implicit none
    save
    type, public, extends(species_c) :: gas_c
        real(kind=8) :: crit_temp !> critical temperature
        real(kind=8) :: crit_press !> critical pressure
        real(kind=8) :: acentric_fact !> acentric factor
    contains
    !> Set
        procedure, public :: set_crit_temp
        procedure, public :: set_crit_press
        procedure, public :: set_acentric_fact
    end type
    
    contains
        subroutine set_crit_temp(this,crit_temp)
            implicit none
            class(gas_c) :: this
            real(kind=8), intent(in) :: crit_temp
            this%crit_temp=crit_temp
        end subroutine
        
        subroutine set_crit_press(this,crit_press)
            implicit none
            class(gas_c) :: this
            real(kind=8), intent(in) :: crit_press
            this%crit_press=crit_press
        end subroutine
        
        subroutine set_acentric_fact(this,acentric_fact)
            implicit none
            class(gas_c) :: this
            real(kind=8), intent(in) :: acentric_fact
            this%acentric_fact=acentric_fact
        end subroutine
        
end module