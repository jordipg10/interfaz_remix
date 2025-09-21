!> Updates indices of equilibrium reactions in reactive zone object
subroutine update_eq_reactions(this,old_eq_reacts_ind)
    use reactive_zone_Lagr_m, only: reactive_zone_c, eq_reaction_c
    implicit none
    class(reactive_zone_c) :: this
    integer(kind=4), intent(in) :: old_eq_reacts_ind(:) !> indices of old equilibrium reactions
    
    integer(kind=4) :: i,j,k
    type(eq_reaction_c), allocatable :: aux(:) !> auxiliary variable equilibrium reactions
    logical :: flag
    integer(kind=4), allocatable :: aux_ind(:)

    aux_ind=this%ind_eq_reacts
    deallocate(this%ind_eq_reacts)
    call this%allocate_ind_eq_reacts()
    j=1
    k=1
    
    do 
        i=1
        flag=.true.
        do
            if (old_eq_reacts_ind(i)==j) then
                flag=.false.
                exit
            else if (i<size(old_eq_reacts_ind)) then
                i=i+1
            else
                exit
            end if
        end do
        if (flag.eqv..true.) then
            this%chem_syst%eq_reacts(this%ind_eq_reacts(k))=aux(j)
            if (k<this%speciation_alg%num_eq_reactions) then
                k=k+1
            end if
        end if
        if (j<size(aux)) then
            j=j+1
        end if
    end do
    
end subroutine update_eq_reactions
