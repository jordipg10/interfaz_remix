!> This subroutine sets the solid stoichiometric matrix of the chemical system from its reactions
subroutine set_stoich_mat_sol(this)
    use chem_system_m, only: chem_system_c
    implicit none
    class(chem_system_c) :: this
    
    integer(kind=4) :: i,j,k,l
    logical :: flag

    if (this%num_reacts>0) then
        if (allocated(this%stoich_mat_sol)) then
            deallocate(this%stoich_mat_sol)
        end if
        allocate(this%stoich_mat_sol(this%num_reacts,this%num_solids))
        this%stoich_mat_sol=0d0 !> initialisation
    !> We set surface complexes
        if (this%cat_exch%num_surf_compl>0) then
            i=1 !> counter total surface complexes
            j=1 !> counter equilibrium reactions reactive zone
            k=1 !> counter species of each equilibrium reaction in reactive zone
            do
                if (this%cat_exch%surf_compl(i)%name==this%eq_reacts(j)%species(k)%name) then
                    this%stoich_mat_sol(j,i)=this%eq_reacts(j)%stoichiometry(k)
                    if (j<this%num_eq_reacts) then
                        j=j+1 
                        k=1
                    else if (i<this%cat_exch%num_surf_compl) then
                        i=i+1
                        j=1
                        k=1
                    else
                        exit
                    end if
                else if (k<this%eq_reacts(j)%num_species) then
                    k=k+1
                else if (j<this%num_eq_reacts) then
                    j=j+1
                    k=1
                else if (i<this%cat_exch%num_surf_compl) then
                    i=i+1
                    j=1
                    k=1
                else
                    exit
                end if
            end do
        end if
    !> We set minerals
        if (this%num_minerals_eq>0) then
            i=1 !> counter total minerals eq
            j=1 !> counter equilibrium reactions reactive zone
            k=1 !> counter species of each equilibrium reaction in chyemical system
            do
                if (this%minerals(this%num_minerals_kin+i)%mineral%name==this%eq_reacts(j)%species(k)%name) then
                    this%stoich_mat_sol(j,this%cat_exch%num_surf_compl+this%num_minerals_kin+i)=this%eq_reacts(j)%stoichiometry(k)
                    if (j<this%num_eq_reacts) then
                        j=j+1 
                        k=1
                    else if (i<this%num_minerals_eq) then
                        i=i+1
                        j=1
                        k=1
                    else
                        exit
                    end if
                else if (k<this%eq_reacts(j)%num_species) then
                    k=k+1
                else if (j<this%num_eq_reacts) then
                    j=j+1
                    k=1
                else if (i<this%num_minerals_eq) then
                    i=i+1
                    j=1
                    k=1
                else
                    exit
                end if
            end do
        end if
        if (this%num_minerals_kin>0) then
            i=1 !> counter total kinetic minerals
            j=1 !> counter mineral kinetic reactions
            k=1 !> counter species of each mineral kinetic reaction
            do
                if (this%minerals(i)%mineral%name==this%min_kin_reacts(j)%species(k)%name) then
                    this%stoich_mat_sol(this%num_eq_reacts+this%num_lin_kin_reacts+j,this%cat_exch%num_surf_compl+i)=&
                    this%min_kin_reacts(j)%stoichiometry(k)
                    if (j<this%num_minerals_kin) then
                        j=j+1
                        k=1
                    else if (i<this%num_minerals_kin) then
                        i=i+1
                        j=1
                        k=1
                    else
                        exit
                    end if
                else if (k<this%min_kin_reacts(j)%num_species) then
                    k=k+1
                else if (j<this%num_minerals_kin) then
                    j=j+1
                    k=1
                else if (i<this%num_minerals_kin) then
                    i=i+1
                    j=1
                    k=1
                else
                    exit
                end if
            end do
        end if
    end if
end subroutine