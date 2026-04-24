!> Gets indices of reaction species in aqueous phase
!> Note: This function uses species_ind which are chemical system indices.
!> It searches for each reaction species in the aqueous phase by checking
!> if the species_ind matches any aqueous species position (1..num_species-wat_flag).
function get_indices_reaction(this,reaction) result(indices)
    use aq_phase_m, only: aq_phase_c
    use reaction_m, only: reaction_c
    implicit none
    class(aq_phase_c), intent(in) :: this !> aqueous chemistry object
    class(reaction_c), intent(in) :: reaction !> reaction object
    integer(kind=4), allocatable :: indices(:) !> indices of reaction species in aqueous phase
    
    integer(kind=4) :: i,j
!> If reaction is mineral dissolution/precipitation or gas in equilibrium, all species are aqueous except one
    if (reaction%react_type==2 .or. reaction%react_type==6) then !> mineral or gas reaction
        allocate(indices(reaction%num_species-1))
    else if (reaction%react_type==3) then !> cation exchange reaction
        allocate(indices(reaction%num_species-2))
    else if (reaction%react_type==4) then !> redox 
        allocate(indices(reaction%num_species-this%wat_flag)) !> (chapuza)
    end if
    i=1 !> counter aqueous species in reaction
    j=1 !> counter species in aqueous phase
    do
        !> species_ind stores indices into the chemical system species array
        !> For aqueous species, this index coincides with the aq_phase index
        if (reaction%species_ind(i)==j) then
            indices(i)=j
            if (i<size(indices)) then
                i=i+1
                j=1
            else
                exit
            end if
        else if (j<this%num_species-this%wat_flag) then
            j=j+1
        else
            error stop "Species is not in aqueous phase"
        end if
    end do
end function