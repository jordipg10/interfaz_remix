subroutine link_target_waters_reactive_zone(this,i,dom_indices,ext_indices)
    use chemistry_m,   only: chemistry_c
    use array_ops_m, only: append_int_1D_array
    implicit none
    
    class(chemistry_c), intent(in) :: this
    integer(kind=4), intent(in) :: i !> reactive zone index
    integer(kind=4), intent(out), allocatable :: dom_indices(:) !> indices of domain target waters associated to i_th reactive zone
    integer(kind=4), intent(out), allocatable :: ext_indices(:) !> indices of external target waters associated to i_th reactive zone
    
    integer(kind=4) :: j,k,l,n,nf_sp_ind
    logical :: flag
    real(kind=8), parameter :: eps=1d-12
    
    allocate(dom_indices(0),ext_indices(0))

    j=1
    k=1
    l=1
    flag=.true.
    !print *, this%dom_tar_wat_indices
    do
        if (this%target_waters(this%dom_tar_wat_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species>0 .and. &
        this%target_waters(this%dom_tar_wat_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==&
        this%reactive_zones(i)%num_non_flowing_species) then
            call this%target_waters(this%dom_tar_wat_indices(j))%solid_chemistry%reactive_zone%is_nf_species_in_react_zone(&
            this%reactive_zones(i)%non_flowing_species(k),flag,nf_sp_ind)
            if (.not. flag) then
                if (j<this%num_target_waters_dom) then
                    j=j+1
                    k=1
                else
                    exit
                end if
            else if (k<this%reactive_zones(i)%num_non_flowing_species) then
                k=k+1
            else if (j<this%num_target_waters_dom) then
                call append_int_1D_array(dom_indices,this%dom_tar_wat_indices(j))
                !if (j==this%ext_wat_indices(l)) then
                !    call append_int_1D_array(dom_indices,j)
                !end if
                j=j+1
            else
                call append_int_1D_array(dom_indices,this%dom_tar_wat_indices(j))
                exit
            end if
        else if (this%target_waters(this%dom_tar_wat_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==0 .and. &
            this%target_waters(this%dom_tar_wat_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==&
            this%reactive_zones(i)%num_non_flowing_species) then
            if (j<this%num_target_waters_dom) then
                call append_int_1D_array(dom_indices,this%dom_tar_wat_indices(j))
                j=j+1
            else
                call append_int_1D_array(dom_indices,this%dom_tar_wat_indices(j))
                exit
            end if
        else if (j<this%num_target_waters_dom) then
            j=j+1
        else
            exit
        end if
    end do
    
    if (this%num_ext_waters>0) then
        j=1
        k=1
        l=1
        flag=.true.
        !print *, this%ext_waters_indices
        do
            
            !print *, this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species
            !print *, this%reactive_zones(i)%num_non_flowing_species
            if (this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species>0 .and. &
            this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==&
            this%reactive_zones(i)%num_non_flowing_species) then
                call this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%is_nf_species_in_react_zone(&
                this%reactive_zones(i)%non_flowing_species(k),flag,nf_sp_ind)
                if (flag.eqv..false.) then
                    if (j<this%num_ext_waters) then
                        j=j+1
                        k=1
                    else
                        exit
                    end if
                else if (k<this%reactive_zones(i)%num_non_flowing_species) then
                    k=k+1
                else if (j<this%num_ext_waters) then
                    call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
                    !if (j==this%ext_wat_indices(l)) then
                    !    call append_int_1D_array(dom_indices,j)
                    !end if
                    j=j+1
                else
                    call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
                    exit
                end if
            else if (this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==0 .and.&
                 this%target_waters(this%ext_waters_indices(j))%solid_chemistry%reactive_zone%num_non_flowing_species==&
                 this%reactive_zones(i)%num_non_flowing_species) then
                if (j<this%num_ext_waters) then
                    call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
                    j=j+1
                else
                    call append_int_1D_array(ext_indices,this%ext_waters_indices(j))
                    exit
                end if
            else if (j<this%num_ext_waters) then
                j=j+1
            else
                exit
            end if
        end do
    end if

end subroutine