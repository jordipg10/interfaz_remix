!> This subroutine finds all waters whose CURRENT state belongs to the i-th reactive zone, i.e.
!> whose solid_chemistry%reactive_zone matches this%reactive_zones(i) (by reactive-zone id, which
!> must be set beforehand): the grouping of the waters by the zone of the cell they sit in. In the
!> Eulerian solvers this is the fixed cell grouping; after a Lagrangian displacement a parcel is
!> grouped with the zone of the cell it has moved INTO (its previous-step state may have been
!> equilibrated in another zone: the within-zone mixing in the CURRENT zone, mix_within_zone_new_rz,
!> brings it to equilibrium with the new zone). The external waters (ext_waters_indices: boundary and
!> recharge) are never grouped in a zone column, even when their solid chemistry carries a domain zone
!> (the external twins of the Dirichlet-node boundary waters of the edge-centred scheme, which
!> set_target_waters copies from the node waters): the mixer receives the same external list and
!> mixes those waters directly, and its partition check requires every stencil water to be in exactly
!> one of {zone column, externals}. Boundary waters without a target solid carry no domain zone
!> (id 0) and would be in no column anyway. Grouping by the previous-step zone
!> (solid_chemistry_old, the companion of mix_within_zone_old_rz) is the variant this replaced.
!> It is the reactive-zone analogue of link_waters_mineral_zone.
subroutine set_wat_indices_zones(this)
    use chemistry_m, only: chemistry_c
    implicit none
    
    class(chemistry_c) :: this !> chemistry object
    !integer(kind=4), intent(in) :: i !> reactive zone index
    !integer(kind=4), intent(out), allocatable :: tar_wat_indices_rz(:) !> global indices of target waters whose previous-step state belongs to the i-th reactive zone
    !integer(kind=4), intent(in), optional :: exit_wat_indices(:) !> global indices of the boundary waters holding parcels that left the domain (Neumann outflow): classified like the target waters
    
    integer(kind=4) :: i,j,n,n_exit
    integer(kind=4), allocatable :: buffer(:)
    logical, allocatable :: ext_mask(:) !> ext_mask(j)=.true. if water j is listed in ext_waters_indices (never grouped in a zone column)
    
    n_exit=0
    !if (present(exit_wat_indices)) n_exit=size(exit_wat_indices)
    !> Single pass into a buffer sized to the number of candidate waters, then trim
    !> once. Avoids the O(n^2) repeated reallocation of appending one index at a time.
    allocate(buffer(this%num_waters))
    !> external waters: mixed directly by the mixer, never grouped in a zone column
    allocate(ext_mask(this%num_waters))
    ext_mask=.false.
    if (allocated(this%ext_waters_indices)) then
        do j=1,size(this%ext_waters_indices)
            if (this%ext_waters_indices(j)>=1 .and. this%ext_waters_indices(j)<=this%num_waters) &
                ext_mask(this%ext_waters_indices(j))=.true.
        end do
    end if
    call this%wat_indices_zones%allocate_array(this%num_reactive_zones)
    do i=1,this%num_reactive_zones
        buffer=0 !> reset buffer for the current reactive zone
        n=0 !> counter for the number of waters in the current reactive zone
        !> waters whose CURRENT reactive zone is this one (match by reactive-zone id, set beforehand)
        do j=1,this%num_waters
            if (ext_mask(j)) cycle
            if (this%waters(j)%solid_chemistry%reactive_zone%id==this%reactive_zones(i)%id) then
                n=n+1
                buffer(n)=j
            end if
        end do
        call this%wat_indices_zones%cols(i)%allocate_vector(n)
        this%wat_indices_zones%cols(i)%col_1=buffer(1:n)
    ! !> exited parcels (outflow boundary waters) whose old reactive zone is this one
    ! do j=1,n_exit
    !     if (this%waters(exit_wat_indices(j))%solid_chemistry_old%reactive_zone%id==this%reactive_zones(i)%id) then
    !         n=n+1
    !         buffer(n)=exit_wat_indices(j)
    !     end if
    ! end do
    !tar_wat_indices_rz=buffer(1:n)
    !this%wat_indices_rz=buffer(1:n)
    end do
end subroutine set_wat_indices_zones