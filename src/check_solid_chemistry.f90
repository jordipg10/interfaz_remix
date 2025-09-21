!> This subroutine checks for zero concentrations in a solid chemistry object
subroutine check_solid_chemistry(this,tolerance,flag,indices)
    use solid_chemistry_m
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(in) :: tolerance !> tolerance for concentrations of solids
    integer(kind=4), intent(out) :: flag !> 1 if no zero concentrations, 0 otherwise
    integer(kind=4), intent(out), allocatable :: indices(:) !> indices of zero concentrations
    
    integer(kind=4) :: num_new_non_flowing_species,i,j,k
    integer(kind=4), allocatable :: old_nf_ind(:),old_solid_ind(:)
    type(species_c), allocatable :: new_non_flowing_species(:)
    type(solid_chemistry_c), allocatable :: new_solid_chems(:)
    real(kind=8), parameter :: epsilon=1d-9 !> arbitrario
    
    flag=1
    do i=1,this%reactive_zone%num_solids
        if (this%concentrations(i)<tolerance) then
            flag=0
            indices=[indices,j] !> indices of zero concentration solids
        else
            continue
        end if
    end do
end subroutine