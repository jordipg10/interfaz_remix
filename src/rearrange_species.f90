subroutine rearrange_species(species,num_prim_species)
    use species_m, only: species_c
    implicit none
    class(species_c) :: species(:)
    integer(kind=4), intent(in), optional :: num_prim_species
    
    integer(kind=4) :: i,j,n
    real(kind=8), parameter :: epsilon=1d-9
    
    n=size(species)

end subroutine