!> Aqueous solution parameters type
!! This type contains parameters related to properties of aqueous solution
module params_aq_sol_m
    implicit none
    save
    type, public :: params_aq_sol_t !> aqueous solution parameters type
        real(kind=8) :: A=0.5092 !> (25ºC)
        real(kind=8) :: B=0.3283 !> (25ºC)
    end type
end module