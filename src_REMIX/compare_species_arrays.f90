!> Checks if two arrays of species are equal or not
subroutine are_species_arrays_equal(species_array_1,species_array_2,flag)
    use species_m
    implicit none
    class(species_c), intent(in) :: species_array_1(:)
    class(species_c), intent(in) :: species_array_2(:)
    logical, intent(out) :: flag !> TRUE if arrays are the same, FALSE otherwise

    integer(kind=4) :: i,n
    
    flag=.true.
    n=size(species_array_1)
    if(size(species_array_1)/=n) error stop "Dimension error in compare_species_arrays"
    
    if (n>0) then
        do i=1,n
            if (species_array_1(i)%name/=species_array_2(i)%name) then
                flag=.false.
                exit
            else
                continue
            end if
        end do
    end if
end subroutine