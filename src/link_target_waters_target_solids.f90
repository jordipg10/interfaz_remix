subroutine link_waters_target_solids(this,tar_sol_indices,wat_indices)
    use chemistry_m,   only: chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use array_ops_m, only: append_int_1D_array
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: tar_sol_indices(:) !> target solid indices
    integer(kind=4), intent(out), allocatable :: wat_indices(:) !> indices of target waters associated to target solids
    
    integer(kind=4) :: i,j,k,n_tw_sol
    logical :: flag
    real(kind=8), parameter :: eps=1d-12

    if (size(tar_sol_indices)==0) then
        allocate(wat_indices(0))
        return
    end if
    
        allocate(wat_indices(size(tar_sol_indices))) !> we assume bijection
    
    
        i=1
        j=1
        do
            if (this%waters(j)%solid_chemistry%reactive_zone%num_non_flow_species==&
                this%target_solids(tar_sol_indices(i))%reactive_zone%num_non_flow_species) then
                if (this%waters(j)%solid_chemistry%reactive_zone%num_non_flow_species>0) then
                    flag = all(this%waters(j)%solid_chemistry%reactive_zone%ind_non_flow_species==&
                        this%target_solids(tar_sol_indices(i))%reactive_zone%ind_non_flow_species)
                else
                    flag = .true.
                end if
                if (flag) then
                    if (allocated(this%waters(j)%solid_chemistry%concentrations) .and. &
                        allocated(this%target_solids(tar_sol_indices(i))%concentrations)) then
                        flag = inf_norm_vec_real(this%waters(j)%solid_chemistry%concentrations-&
                            this%target_solids(tar_sol_indices(i))%concentrations)<eps
                    end if
                end if
                if (flag) then
                    wat_indices(i)=j
                    if (j<this%num_waters) then
                        j=j+1
                    end if
                    if (i<size(tar_sol_indices)) then
                        i=i+1
                    else
                        exit
                    end if
                else if (j<this%num_waters) then
                    j=j+1
                else
                    exit
                end if
            else if (j<this%num_waters) then
                j=j+1
            !else if (i<size(tar_sol_indices)) then
            !>    i=i+1
            !>    !i=1
            else
                error stop "Error in link_waters_target_solids"
            end if
        end do
end subroutine