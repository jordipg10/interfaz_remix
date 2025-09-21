!> Mobile zone module
module mob_zone_m
    implicit none
    save
    type, public :: mob_zone_c
        !real(kind=8), allocatable :: conc(:)
        real(kind=8) :: mob_por !> mobile porosity
    contains
        !procedure, public :: set_conc_init
    end type
    
    contains
        !subroutine set_conc_init(this,conc_init)
        !>    implicit none
        !>    class(mob_zone_c) :: this
        !>    real(kind=8), intent(in) :: conc_init(:)
        !>    this%conc=conc_init
        !end subroutine
end module