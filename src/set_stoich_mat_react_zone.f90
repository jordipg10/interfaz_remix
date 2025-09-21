!> This subroutine sets the stoichiometric matrix of the reactive zone from its reactions
subroutine set_stoich_mat_react_zone(this)
    use reactive_zone_Lagr_m, only: reactive_zone_c
    implicit none
    class(reactive_zone_c) :: this
    
    integer(kind=4) :: i,j,k,l,n_k
    logical :: flag
            
    i=1 !> counter total species chemical system
    j=1 !> counter equilibrium reactions reactive zone
    k=1 !> counter species of each equilibrium reaction in reactive zone
    l=1 !> counter total species reactive zone
    if (this%speciation_alg%num_eq_reactions>0) then
        if (allocated(this%stoich_mat)) then
            deallocate(this%stoich_mat)
        end if
        allocate(this%stoich_mat(this%speciation_alg%num_eq_reactions,this%speciation_alg%num_species))
        this%stoich_mat=0d0 !> initialisation
        do
            if (this%chem_syst%species(i)%name==this%chem_syst%eq_reacts(this%ind_eq_reacts(j))%species(k)%name) then
                flag=.true. !> species is involved in reactive zone
                this%stoich_mat(j,l)=this%chem_syst%eq_reacts(this%ind_eq_reacts(j))%stoichiometry(k)
                if (j<this%speciation_alg%num_eq_reactions) then
                    j=j+1 
                    k=1
                else if (l<this%speciation_alg%num_species) then
                    i=i+1
                    l=l+1
                    j=1
                    k=1
                    flag=.false.
                else
                    exit
                end if
            else if (k<this%chem_syst%eq_reacts(this%ind_eq_reacts(j))%num_species) then
                k=k+1
            else if (j<this%speciation_alg%num_eq_reactions) then
                j=j+1
                k=1
            else if (i<this%chem_syst%num_species .and. this%speciation_alg%flag_comp .eqv. .true.) then
                if ((flag .eqv. .true.) .or. (i<=this%speciation_alg%num_species-this%num_minerals-&
                    this%gas_phase%num_gases_eq)) then
                    if (l<this%speciation_alg%num_species) then
                        l=l+1
                    else
                        exit
                    end if
                end if
                i=i+1
                j=1
                k=1
                flag=.false.
            else if (i<this%chem_syst%num_species .and. this%speciation_alg%flag_comp .eqv. .false.) then
                if ((flag .eqv. .true.) .or. (i<=this%speciation_alg%num_aq_prim_species)) then
                    if (l<this%speciation_alg%num_species) then
                        l=l+1
                    else
                        exit
                    end if
                end if
                i=i+1
                j=1
                k=1
                flag=.false.
            else
                exit
            end if
        end do
    end if
end subroutine