!> Appends reaction object to array of reaction objects
subroutine append_reaction_to_array(this,reaction)
    use chem_system_m
    implicit none
    class(chem_system_c) :: this
    class(reaction_c), intent(in) :: reaction
    
    integer(kind=4) :: i
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    type(redox_kin_c), allocatable :: redox_kin_reacts(:)
    type(kin_reaction_poly_c), allocatable :: kin_react_ptrs(:)
    
    select type (reaction)
    type is (eq_reaction_c)
        eq_reacts=this%eq_reacts
        if (allocated(this%eq_reacts)) then
            deallocate(this%eq_reacts)
        end if
        call this%allocate_eq_reacts(size(eq_reacts)+1)
        do i=1,this%num_eq_reacts-1
            this%eq_reacts(i)=eq_reacts(i)
        end do
        this%eq_reacts(this%num_eq_reacts)=reaction
    type is (redox_kin_c)
        redox_kin_reacts=this%redox_kin_reacts
        if (allocated(this%redox_kin_reacts)) then
            deallocate(this%redox_kin_reacts)
        end if
        call this%allocate_redox_kin_reacts(size(redox_kin_reacts)+1)
        do i=1,this%num_redox_kin_reacts-1
            this%redox_kin_reacts(i)=redox_kin_reacts(i)
        end do
        this%redox_kin_reacts(this%num_redox_kin_reacts)=reaction
    class default
        error stop "Allocate subroutine not implemented for this reaction subclass yet"
    end select
end subroutine append_reaction_to_array