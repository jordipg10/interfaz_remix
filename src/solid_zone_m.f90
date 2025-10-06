module solid_zone_m
    use solid_m
    implicit none
    save
    type, public :: solid_zone_c
        type(solid_c), allocatable :: solids(:)
        integer(kind=4) :: num_solids
    contains
    end type
end module