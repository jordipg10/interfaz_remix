subroutine link_target_waters_target_gases(this,tar_gas_indices,tar_wat_indices)
    use chemistry_m,  only: chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use array_ops_m, only: append_int_1D_array
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: tar_gas_indices(:) !> target gas indices
    integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) !> indices of target waters associated to target gases
    
    integer(kind=4) :: i,j,k,n_tw_sol
    logical :: flag
    real(kind=8), parameter :: eps=1d-12

    if (size(tar_gas_indices)>0) then
    
        allocate(tar_wat_indices(size(tar_gas_indices))) !> we assume bijection
    
    
        i=1
        j=1
        do
            if (inf_norm_vec_real(this%target_waters(j)%gas_chemistry%concentrations-this%target_gases(&
            tar_gas_indices(i))%concentrations)<eps) then
                tar_wat_indices(i)=j
                if (j<this%num_target_waters) then
                    j=j+1
                end if
                if (i<size(tar_gas_indices)) then
                    i=i+1
                else
                    exit
                end if
            else if (j<this%num_target_waters) then
                j=j+1
            else
                error stop "Error in link_target_waters_target_gases"
            end if
        end do
    end if
end subroutine