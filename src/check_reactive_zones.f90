!!> This subroutine checks if reactive zones have changed
!subroutine check_reactive_zones(this)
!>    use chemistry_Lagr_m
!>    implicit none
!>    class(chemistry_c) :: this
!>    
!>    integer(kind=4) :: num_new_react_zones,i,j
!>    integer(kind=4), allocatable :: old_nf_ind(:)
!>    type(species_c), allocatable :: new_non_flowing_species(:)
!>    type(reactive_zone_c) :: new_react_zone
!>    real(kind=8), parameter :: epsilon=1d-9 !> arbitrary
!>    
!>    do i=1,this%num_target_solids
!>        do j=1,this%target_solids(i)%reactive_zone%num_non_flowing_species
!>            if (this%target_solids(i)%concentrations(j)<epsilon) then
!>                num_new_react_zones=num_new_react_zones+1
!>                !> Autentica chapuza
!>                old_nf_ind=[old_nf_ind,j] !> indices of non flowing species that not belong anymore to reactive zone
!>            else
!>                !> Autentica chapuza
!>                new_non_flowing_species=[new_non_flowing_species,this%target_solids(i)%reactive_zone%non_flowing_species(j)]
!>            end if
!>        end do
!>        new_react_zone=this%target_solids(i)%reactive_zone !> we initialise the new reactive zone
!>        !call new_react_zone%set_non_flowing_species(new_non_flowing_species) !> we set the new non flowing species
!>        call new_react_zone%update_reactive_zone(old_nf_ind,new_react_zone) !> we update the attributes of the new reactive zone
!>        !> Autentica chapuza
!>        this%reactive_zones=[this%reactive_zones,new_react_zone]
!>        this%num_reactive_zones=this%num_reactive_zones+1
!>        call this%target_solids(i)%set_reactive_zone(this%reactive_zones(this%num_reactive_zones))
!>    end do
!end subroutine