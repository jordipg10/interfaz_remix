subroutine link_waters_mineral_zone(this,i,wat_indices)
    use chemistry_m, only: chemistry_c
    use mineral_zone_m, only: mineral_zone_c
    use array_ops_m, only: append_int_1D_array
    implicit none
    
    class(chemistry_c), intent(in) :: this
    integer(kind=4), intent(in) :: i !> mineral zone index
    integer(kind=4), intent(out), allocatable :: wat_indices(:) !> indices of domain target waters associated to i_th mineral zone
    
    integer(kind=4) :: j,k,l,n,min_ind
    logical :: flag
    real(kind=8), parameter :: eps=1d-16
    
    allocate(wat_indices(0))

    j=1
    k=1
    l=1
    flag=.true.
    do
        if (this%waters(j)%solid_chemistry%mineral_zone%num_minerals>0 .and. &
            this%waters(j)%solid_chemistry%mineral_zone%num_minerals==&
            this%mineral_zones(i)%num_minerals) then
            call this%waters(j)%solid_chemistry%mineral_zone%is_mineral_in_min_zone(&
                this%mineral_zones(i)%chem_syst%minerals(this%mineral_zones(i)%ind_min_chem_syst(k)),flag,min_ind)
            if (flag .eqv. .false.) then
                if (j<this%num_waters) then
                    j=j+1
                    k=1
                else
                    exit
                end if
            else if (k<this%mineral_zones(i)%num_minerals) then
                k=k+1
            else if (j<this%num_waters) then
                call append_int_1D_array(wat_indices,j) !> chapuza
                !if (j==this%ext_wat_indices(l)) then
                !    call append_int_1D_array(wat_indices,j)
                !end if
                j=j+1
            else
                call append_int_1D_array(wat_indices,j) !> chapuza
                exit
            end if
        else if (this%waters(j)%solid_chemistry%mineral_zone%num_minerals==0 .and. &
            this%waters(j)%solid_chemistry%mineral_zone%num_minerals==&
            this%mineral_zones(i)%num_minerals) then
            if (j<this%num_waters) then
                call append_int_1D_array(wat_indices,j) !> chapuza
                j=j+1
            else
                call append_int_1D_array(wat_indices,j) !> chapuza
                exit
            end if
        else if (j<this%num_waters) then
            j=j+1
        else
            exit
        end if
    end do
    
    !if (this%num_ext_waters>0) then
    !    j=1
    !    k=1
    !    l=1
    !    flag=.true.
    !    do
    !        if (this%waters(this%ext_waters_indices(j))%solid_chemistry%mineral_zone%num_minerals>0 .and. &
    !        this%waters(this%ext_waters_indices(j))%solid_chemistry%mineral_zone%num_minerals==&
    !        this%mineral_zones(i)%num_minerals) then
    !            call this%waters(this%ext_waters_indices(j))%solid_chemistry%mineral_zone%is_mineral_in_min_zone(&
    !                this%mineral_zones(i)%chem_syst%minerals(this%mineral_zones(i)%ind_min_chem_syst(k)),flag,min_ind)
    !            if (flag .eqv. .false.) then
    !                if (j<this%num_ext_waters) then
    !                    j=j+1
    !                    k=1
    !                else
    !                    exit
    !                end if
    !            else if (k<this%mineral_zones(i)%num_minerals) then
    !                k=k+1
    !            else if (j<this%num_ext_waters) then
    !                call append_int_1D_array(ext_indices,this%ext_waters_indices(j)) !> chapuza
    !                !if (j==this%ext_wat_indices(l)) then
    !                !    call append_int_1D_array(wat_indices,j)
    !                !end if
    !                j=j+1
    !            else
    !                call append_int_1D_array(ext_indices,this%ext_waters_indices(j))    !> chapuza
    !                exit
    !            end if
    !        else if (this%waters(this%ext_waters_indices(j))%solid_chemistry%mineral_zone%num_minerals==0 .and. &
    !            this%waters(this%ext_waters_indices(j))%solid_chemistry%mineral_zone%num_minerals==&
    !            this%mineral_zones(i)%num_minerals) then
    !            if (j<this%num_ext_waters) then
    !                call append_int_1D_array(ext_indices,this%ext_waters_indices(j)) !> chapuza
    !                j=j+1
    !            else
    !                call append_int_1D_array(ext_indices,this%ext_waters_indices(j)) !> chapuza
    !                exit
    !            end if
    !        else if (j<this%num_ext_waters) then
    !            j=j+1
    !        else
    !            exit
    !        end if
    !    end do
    !end if

end subroutine