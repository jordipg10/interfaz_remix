!> This subroutine checks if two reactive zones have the same non-flowing species
subroutine compare_react_zones(react_zone_1,react_zone_2,flag)
    use reactive_zone_Lagr_m, only: reactive_zone_c
    implicit none
    class(reactive_zone_c), intent(in) :: react_zone_1 !> first reactive zone
    class(reactive_zone_c), intent(in) :: react_zone_2 !> second reactive zone
    logical, intent(out) :: flag !> TRUE if same non flowing species, FALSE otherwise

    integer(kind=4) :: i,n,nf_species_ind
    logical :: sp_flag
    
    flag=.true.
    if (react_zone_1%num_non_flowing_species/=react_zone_2%num_non_flowing_species) then
        flag=.false.
    else
        do i=1,react_zone_1%num_non_flowing_species
            call react_zone_2%is_nf_species_in_react_zone(react_zone_1%non_flowing_species(i),sp_flag,nf_species_ind)
            if (sp_flag .eqv. .false.) then
                flag=.false.
                exit
            end if
        end do
    end if
end subroutine