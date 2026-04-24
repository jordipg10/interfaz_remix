subroutine link_target_solids_reactive_zone(this,i,tar_sol_indices)
    use chemistry_m, only : chemistry_c
    use array_ops_m, only : append_int_1D_array
    implicit none
    
    class(chemistry_c), intent(in) :: this !> chemistry object
    integer(kind=4), intent(in) :: i !> reactive zone index
    integer(kind=4), intent(out), allocatable :: tar_sol_indices(:) !> indices of domain target waters associated to i_th reactive zone
    !integer(kind=4), intent(out), allocatable :: ext_indices(:) !> indices of external target waters associated to i_th reactive zone
    
    integer(kind=4) :: j,k,l,n,nf_sp_ind
    logical :: flag
    real(kind=8), parameter :: eps=1d-12
    
    allocate(tar_sol_indices(0))

    j=1
    k=1
    l=1
    flag=.true.
    do
        if (this%target_solids(j)%reactive_zone%num_non_flow_species>0 .and. &
        this%target_solids(j)%reactive_zone%num_non_flow_species==&
        this%reactive_zones(i)%num_non_flow_species) then
            call this%target_solids(j)%reactive_zone%is_nf_species_in_react_zone(&
            this%reactive_zones(i)%ind_non_flow_species(k),flag,nf_sp_ind)
            if (flag.eqv..false.) then
                if (j<this%num_target_solids) then
                    j=j+1
                    k=1
                else
                    exit
                end if
            else if (k<this%reactive_zones(i)%num_non_flow_species) then
                k=k+1
            else if (j<this%num_target_solids) then
                call append_int_1D_array(tar_sol_indices,j)
                !if (j==this%ext_wat_indices(l)) then
                !    call append_int_1D_array(dom_indices,j)
                !end if
                j=j+1
            else
                call append_int_1D_array(tar_sol_indices,j)
                exit
            end if
        else if (this%target_solids(j)%reactive_zone%num_non_flow_species==0 .and. &
            this%target_solids(j)%reactive_zone%num_non_flow_species==&
            this%reactive_zones(i)%num_non_flow_species) then
            if (j<this%num_target_solids) then
                call append_int_1D_array(tar_sol_indices,j)
                j=j+1
            else
                call append_int_1D_array(tar_sol_indices,j)
                exit
            end if
        else if (j<this%num_target_solids) then
            j=j+1
        else
            exit
        end if
    end do
    !j=1
    !k=1
    !l=1
    !flag=.true.
    !do
    !    if (this%target_solids(this%ext_waters_indices(j))%reactive_zone%num_non_flow_species>0 .and. &
    !    this%target_solids(this%ext_waters_indices(j))%reactive_zone%num_non_flow_species==&
    !    this%reactive_zones(i)%num_non_flow_species) then
    !        call this%target_solids(this%ext_waters_indices(j))%reactive_zone%is_nf_species_in_react_zone(&
    !        this%reactive_zones(i)%ind_non_flow_species(k),flag,nf_sp_ind)
    !        if (flag.eqv..false.) then
    !            if (j<this%num_ext_waters) then
    !                j=j+1
    !                k=1
    !            else
    !                exit
    !            end if
    !        else if (k<this%reactive_zones(i)%num_non_flow_species) then
    !            k=k+1
    !        else if (j<this%num_ext_waters) then
    !            call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
    !            !if (j==this%ext_wat_indices(l)) then
    !            !    call append_int_1D_array(dom_indices,j)
    !            !end if
    !            j=j+1
    !        else
    !            call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
    !            exit
    !        end if
    !    else if (this%target_solids(this%ext_waters_indices(j))%reactive_zone%num_non_flow_species==0 .and. &
    !        this%target_solids(this%ext_waters_indices(j))%reactive_zone%num_non_flow_species==&
    !        this%reactive_zones(i)%num_non_flow_species) then
    !        if (j<this%num_ext_waters) then
    !            call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
    !            j=j+1
    !        else
    !            call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
    !            exit
    !        end if
    !    else if (j<this%num_ext_waters) then
    !        j=j+1
    !    else
    !        exit
    !    end if
    !end do

end subroutine