subroutine link_target_solids_mineral_zone(this,i,tar_sol_indices)
    use chemistry_m, only: chemistry_c
    use mineral_zone_m, only: mineral_zone_c
    use array_ops_m, only: append_int_1D_array
    implicit none
    
    class(chemistry_c), intent(in) :: this
    integer(kind=4), intent(in) :: i !> mineral zone index
    integer(kind=4), intent(out), allocatable :: tar_sol_indices(:) !> indices of target solids associated to i_th mineral zone
    
    integer(kind=4) :: j,k,min_ind
    logical :: flag
    
    allocate(tar_sol_indices(0))

    j=1
    k=1
    flag=.true.
    do
        if (this%target_solids(j)%mineral_zone%num_minerals>0 .and. &
            this%target_solids(j)%mineral_zone%num_minerals==&
            this%mineral_zones(i)%num_minerals) then
            call this%target_solids(j)%mineral_zone%is_mineral_in_min_zone(&
                this%mineral_zones(i)%chem_syst%minerals(this%mineral_zones(i)%ind_min_chem_syst(k)),flag,min_ind)
            if (flag .eqv. .false.) then
                if (j<this%num_target_solids) then
                    j=j+1
                    k=1
                else
                    exit
                end if
            else if (k<this%mineral_zones(i)%num_minerals) then
                k=k+1
            else if (j<this%num_target_solids) then
                call append_int_1D_array(tar_sol_indices,j)
                j=j+1
                k=1
            else
                call append_int_1D_array(tar_sol_indices,j)
                exit
            end if
        else if (this%target_solids(j)%mineral_zone%num_minerals==0 .and. &
            this%mineral_zones(i)%num_minerals==0) then
            if (j<this%num_target_solids) then
                call append_int_1D_array(tar_sol_indices,j)
                j=j+1
            else
                call append_int_1D_array(tar_sol_indices,j)
                exit
            end if
        else if (j<this%num_target_solids) then
            j=j+1
            k=1
        else
            exit
        end if
    end do
end subroutine
