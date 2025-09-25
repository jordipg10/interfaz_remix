module conc_m
    implicit none
    type, public :: conc_c !> Concentration class for diffusion or transport problems
        real(kind=8), allocatable :: conc(:) !> concentration (c)
        real(kind=8), allocatable :: conc_ext(:) !> external concentration (c_e)
        real(kind=8), allocatable :: conc_bd(:) !> boundary concentration (c_b)
        real(kind=8), allocatable :: conc_r(:) !> recharge concentration (c_r)
        integer(kind=4), allocatable :: conc_r_flag(:)      !> 1 if r>0
                                                            !> 0 otherwise
    contains
        !procedure :: initialise => initialise_conc
        procedure :: set_conc
        procedure :: set_conc_ext
        procedure :: set_conc_bd
        procedure :: set_conc_r
        !procedure :: get_conc => get_conc
    end type conc_c
    
    contains
    
    subroutine set_conc(this, conc)
        class(conc_c) :: this
        real(kind=8), intent(in) :: conc(:)
        this%conc = conc
    end subroutine set_conc
    
    subroutine set_conc_ext(this, conc_ext)
        class(conc_c) :: this
        real(kind=8), intent(in) :: conc_ext(:)
        this%conc_ext = conc_ext
    end subroutine set_conc_ext
    
    subroutine set_conc_bd(this, conc_bd)
        class(conc_c) :: this
        real(kind=8), intent(in) :: conc_bd(:)
        this%conc_bd = conc_bd
    end subroutine set_conc_bd
    
    subroutine set_conc_r(this, conc_r, conc_r_flag)
        class(conc_c) :: this
        real(kind=8), intent(in) :: conc_r(:)
        integer(kind=4), intent(in) :: conc_r_flag(:)
        this%conc_r = conc_r
        this%conc_r_flag = conc_r_flag
    end subroutine set_conc_r
end module conc_m