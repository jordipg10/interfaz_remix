subroutine link_target_waters_old_reaction_zones(this,num_rz,tar_wat_indices_rz)
    use chemistry_m, only: chemistry_c
    use arrays_m, only: int_array_c
    use aqueous_chemistry_m, only: compare_old_reaction_zones
    implicit none
    
    class(chemistry_c), intent(in) :: this !> chemistry object
    integer(kind=4), intent(out) :: num_rz !> number of old reaction zones
    class(int_array_c), intent(out) :: tar_wat_indices_rz !> per old reaction zone (column): global indices of target waters sharing the same reaction set
    
    integer(kind=4) :: i,j,wat_idx
    integer(kind=4), allocatable :: rep(:) !> representative (first) target-water index of each old reaction zone
    integer(kind=4), allocatable :: zone_of(:) !> old reaction-zone index assigned to each target water
    integer(kind=4), allocatable :: count_rz(:) !> number of target waters in each old reaction zone
    integer(kind=4), allocatable :: fill(:) !> running fill counter per old reaction zone
    logical :: flag
    
    !> First pass: assign each target water to an old reaction zone, opening a new zone
    !> whenever the water does not share its full reaction set (equilibrium + kinetic) with
    !> any existing zone representative. At most num_target_waters distinct zones can appear.
    allocate(zone_of(this%num_target_waters))
    allocate(rep(this%num_target_waters))
    allocate(count_rz(this%num_target_waters))
    count_rz=0
    num_rz=0
    do j=1,this%num_target_waters
        wat_idx=this%tar_wat_indices(j)
        zone_of(j)=0
        do i=1,num_rz
            call compare_old_reaction_zones(this%waters(wat_idx),this%waters(rep(i)),flag)
            if (flag) then
                zone_of(j)=i
                exit
            end if
        end do
        if (zone_of(j)==0) then !> water opens a new old reaction zone
            num_rz=num_rz+1
            rep(num_rz)=wat_idx
            zone_of(j)=num_rz
        end if
        count_rz(zone_of(j))=count_rz(zone_of(j))+1
    end do
    
    !> Second pass: one column per old reaction zone, each holding the global indices of
    !> the target waters that belong to that zone.
    call tar_wat_indices_rz%allocate_array(num_rz)
    do i=1,num_rz
        call tar_wat_indices_rz%cols(i)%allocate_vector(count_rz(i))
    end do
    allocate(fill(num_rz))
    fill=0
    do j=1,this%num_target_waters
        i=zone_of(j)
        fill(i)=fill(i)+1
        tar_wat_indices_rz%cols(i)%col_1(fill(i))=this%tar_wat_indices(j)
    end do

end subroutine