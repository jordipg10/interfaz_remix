!> This subroutine checks if two reactive zones have the same non-flowing species.
!> Non-flowing species are stored in chemical-system order, so two equal zones have
!> identical ind_non_flow_species arrays; a direct element-wise comparison suffices
!> (O(n), vectorisable) instead of a per-species search (O(n^2)).
subroutine compare_react_zones(react_zone_1,react_zone_2,flag)
    use reactive_zone_m, only: reactive_zone_c
    implicit none
    class(reactive_zone_c), intent(in) :: react_zone_1 !> first reactive zone
    class(reactive_zone_c), intent(in) :: react_zone_2 !> second reactive zone
    logical, intent(out) :: flag !> TRUE if same non flowing species, FALSE otherwise

    integer(kind=4) :: n

    n=react_zone_1%num_non_flow_species
    if (n/=react_zone_2%num_non_flow_species) then
        flag=.false.
    else if (n==0) then
        flag=.true.
    else
        flag=all(react_zone_1%ind_non_flow_species(1:n)==react_zone_2%ind_non_flow_species(1:n))
    end if
end subroutine