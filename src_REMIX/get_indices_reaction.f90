!> Gets indices of reaction species in aqueous phase
function get_indices_reaction(this,reaction) result(indices)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, reaction_c
    implicit none
    class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry object
    class(reaction_c), intent(in) :: reaction !> reaction object
    integer(kind=4), allocatable :: indices(:) !> indices of reaction species in aqueous phase
    
    integer(kind=4) :: i,j
!> If reaction is mineral dissolution/precipitation or gas in equilibrium, all species are aqueous except one
    if (reaction%react_type==2 .or. reaction%react_type==6) then !> mineral or gas reaction
        allocate(indices(reaction%num_species-1))
    else if (reaction%react_type==3) then !> cation exchange reaction
        allocate(indices(reaction%num_species-2))
    else if (reaction%react_type==4) then !> redox 
        allocate(indices(reaction%num_species-this%aq_phase%wat_flag)) !> (chapuza)
    end if
    i=1 !> counter aqueous species in reaction
    j=1 !> counter species in aqueous phase
    do
        if (reaction%species(i)%name==this%aq_phase%aq_species(j)%name) then
            indices(i)=j
            !print *, reaction%species(i)%name
            if (i<size(indices)) then
                i=i+1
                j=1
            else
                exit
            end if
        else if (j<this%aq_phase%num_species-this%aq_phase%wat_flag) then
            j=j+1
        !else if (i<size(indices)) then
        !    i=i+1
        !    j=1
        else
            error stop "Species is not in aqueous phase"
        end if
    end do
end function