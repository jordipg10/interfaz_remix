!> This subroutine sets the gas stoichiometric matrix of the reactive zone from its reactions
subroutine set_stoich_mat_gas(this)
    use chem_system_m, only: chem_system_c
    implicit none
    class(chem_system_c) :: this
    
    integer(kind=4) :: i,j,k,l
    logical :: flag

!> We set gases
    if (this%gas_phase%num_species>0) then
        if (allocated(this%stoich_mat_gas)) then
            deallocate(this%stoich_mat_gas)
        end if
        allocate(this%stoich_mat_gas(this%num_reacts,this%gas_phase%num_species))
        this%stoich_mat_gas=0d0 !> initialisation
            i=1 !> counter total gases
            j=1 !> counter equilibrium reactions reactive zone
            k=1 !> counter species of each equilibrium reaction in reactive zone
            do
                if (this%gas_phase%gases(i)%name==this%eq_reacts(j)%species(k)%name) then
                    this%stoich_mat_gas(j,i)=this%eq_reacts(j)%stoichiometry(k)
                    if (j<this%num_eq_reacts) then
                        j=j+1 
                        k=1
                    else if (i<this%gas_phase%num_species) then
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
                else if (i<this%gas_phase%num_species) then
                    i=i+1
                    j=1
                    k=1
                else
                    exit
                end if
            end do
            i=1 !> counter total gases
            j=1 !> counter redox reactions reactive zone
            k=1 !> counter species of each redox reaction in reactive zone
            ! CORREGIR EL BUCLE DE ABAJO
            !do
            !    if (this%gas_phase%gases(i)%name==this%redox_kin_reacts(j)%species(k)%name) then
            !        this%stoich_mat_gas(this%num_eq_reacts+j,i)=this%redox_kin_reacts(j)%stoichiometry(k)
            !        if (j<this%num_redox_kin_reacts) then
            !            j=j+1 
            !            k=1
            !        else if (i<this%gas_phase%num_species) then
            !            i=i+1
            !            j=1
            !            k=1
            !        else
            !            exit
            !        end if
            !    else if (k<this%redox_kin_reacts(j)%num_species) then
            !        k=k+1
            !    else if (j<this%num_redox_kin_reacts) then
            !        j=j+1
            !        k=1
            !    else if (i<this%gas_phase%num_species) then
            !        i=i+1
            !        j=1
            !        k=1
            !    else
            !        exit
            !    end if
            !end do
    end if
end subroutine