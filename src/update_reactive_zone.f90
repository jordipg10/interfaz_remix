!> Creates new reactive zone (UNFINISHED)
subroutine update_reactive_zone(this,old_nf_ind)
    use reactive_zone_Lagr_m, only: reactive_zone_c
    implicit none
    class(reactive_zone_c) :: this
    integer(kind=4), intent(in) :: old_nf_ind(:) !> indices non flowing species not present anymore in reactive zone
    
    integer(kind=4) :: i,j,n_old
    integer(kind=4), allocatable :: old_eq_reacts_ind(:),old_min_ind(:)
end subroutine update_reactive_zone