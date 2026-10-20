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
    !> [WAT-ZONES] debug trace (stdout): zone ids, per-zone member lists, waters in no column, parcels
    !> whose old and current zones differ, and a classification summary
    logical, parameter :: dbg=.true. !> set .false. to silence the [WAT-ZONES] prints
    integer(kind=4), parameter :: dbg_max_calls=5 !> full per-zone/per-water listing only for the first calls; summary and warnings afterwards
    integer(kind=4), save :: n_calls=0 !> calls so far (re-run after each Lagrangian displacement when the relinking is enabled)
    integer(kind=4) :: n_tot !> number of waters classified into some zone column
    integer(kind=4) :: n_orphan !> number of waters whose current zone id matches no reactive zone (no column: boundary waters)
    integer(kind=4) :: n_moved !> number of waters whose previous-step and current reactive zones differ
    integer(kind=4) :: old_id,cur_id !> reactive zone ids of a water's previous-step and current solid chemistry
    logical :: is_ext !> .true. if a water is listed among the external (boundary/recharge) waters
    
    n_calls=n_calls+1
    if (dbg) then
        if (n_calls<=dbg_max_calls) then
            print '(A,I0,A,I0,A,I0,A,I0)', "[WAT-ZONES] set_wat_indices_zones call #", n_calls, &
                ": num_waters=", this%num_waters, " num_target_waters=", this%num_target_waters, &
                " num_reactive_zones=", this%num_reactive_zones
            print '(A,*(1X,I0))', "[WAT-ZONES]   reactive zone ids (by position) =", &
                (this%reactive_zones(i)%id, i=1,this%num_reactive_zones)
            !> external (boundary/recharge) waters of the chemistry object: mixed directly by the mixer,
            !> never grouped in a zone column (the Lagrangian solvers extend this list themselves);
            !> listed as water index(water id)
            if (allocated(this%ext_waters_indices)) then
                print '(A,I0,A)', "[WAT-ZONES]   ext_waters_indices (", size(this%ext_waters_indices), "), index(id):"
                if (size(this%ext_waters_indices)>0) print '(("[WAT-ZONES]     ",12(1X,I0,"(",I0,")")))', &
                    (this%ext_waters_indices(j), this%waters(this%ext_waters_indices(j))%id, j=1,size(this%ext_waters_indices))
            else
                print '(A)', "[WAT-ZONES]   ext_waters_indices not allocated"
            end if
        end if
        !> genuine anomalies: printed on every call
        do i=1,this%num_reactive_zones
            if (this%reactive_zones(i)%id/=i) print '(A,I0,A,I0,A)', "[WAT-ZONES]   WARNING: reactive_zones(", i, &
                ")%id=", this%reactive_zones(i)%id, " /= position (cross_zone_mix_core indexes wat_indices_zones%cols by zone id)"
        end do
        !> pointer checks before the classification loop dereferences them
        do j=1,this%num_waters
            if (.not. associated(this%waters(j)%solid_chemistry)) then
                print '(A,I0,A)', "[WAT-ZONES]   WARNING: water ", j, ": solid_chemistry not associated (classification below will fail)"
            else if (.not. associated(this%waters(j)%solid_chemistry%reactive_zone)) then
                print '(A,I0,A)', "[WAT-ZONES]   WARNING: water ", j, &
                    ": solid_chemistry%reactive_zone not associated (classification below will fail)"
            end if
        end do
    end if
    n_tot=0
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
        n_tot=n_tot+n
        if (dbg .and. n_calls<=dbg_max_calls) then
            print '(A,I0,A,I0,A,I0,A)', "[WAT-ZONES]   zone ", i, " (id=", this%reactive_zones(i)%id, "): n_waters=", n, &
                ", members index(id):"
            !> member water index(water id), 12 per line (format reversion repeats the whole group, prefix included)
            if (n>0) print '(("[WAT-ZONES]     ",12(1X,I0,"(",I0,")")))', (buffer(j), this%waters(buffer(j))%id, j=1,n)
        end if
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
    if (dbg) then
        !> waters in no column (CURRENT zone id matching no reactive zone: the boundary waters) and
        !> waters whose previous-step and current zones differ (parcels displaced across an interface)
        n_orphan=0
        n_moved=0
        do j=1,this%num_waters
            cur_id=this%waters(j)%solid_chemistry%reactive_zone%id
            if (.not. any(this%reactive_zones(1:this%num_reactive_zones)%id==cur_id)) then
                n_orphan=n_orphan+1
                is_ext=.false.
                if (allocated(this%ext_waters_indices)) is_ext=any(this%ext_waters_indices==j)
                !> a water in no column must be mixed as an external water by the mixer
                if (n_calls<=dbg_max_calls .or. .not. is_ext) then
                    if (is_ext) then
                        print '(A,I0,A,I0,A)', "[WAT-ZONES]   water ", j, " in NO zone column: current zone id=", cur_id, &
                            " (listed in chemistry%ext_waters_indices)"
                    else
                        print '(A,I0,A,I0,A)', "[WAT-ZONES]   WARNING: water ", j, " in NO zone column: current zone id=", &
                            cur_id, " and NOT in chemistry%ext_waters_indices: the solver must list it among the mixer's" // &
                            " external waters or the partition check fails for stencils containing it"
                    end if
                end if
            end if
            old_id=-1
            if (associated(this%waters(j)%solid_chemistry_old)) then
                if (associated(this%waters(j)%solid_chemistry_old%reactive_zone)) &
                    old_id=this%waters(j)%solid_chemistry_old%reactive_zone%id
            end if
            if (cur_id/=old_id) then
                n_moved=n_moved+1
                if (n_calls<=dbg_max_calls) print '(A,I0,A,I0,A,I0,A)', "[WAT-ZONES]   water ", j, ": previous-step zone id=", &
                    old_id, " current zone id=", cur_id, " (displaced across zones; -1 = previous-step zone not associated)"
            end if
        end do
        if (n_calls<=dbg_max_calls) print '(A,I0,A,I0,A,I0,A,I0,A,I0)', "[WAT-ZONES]   classified=", n_tot, &
            " of num_waters=", this%num_waters, "  externals (never in a column)=", count(ext_mask), &
            "  no-column=", n_orphan, "  previous/current zone differ=", n_moved
    end if
end subroutine set_wat_indices_zones