! !> Updates mineral zone
! subroutine update_mineral_zone(this,old_min_ind)
!     use mineral_zone_m, only: mineral_zone_c, mineral_c
!     class(mineral_zone_c) :: this
!     integer(kind=4), intent(in) :: old_min_ind(:) !> indices of old mineral zone
    
!     integer(kind=4) :: i,j,k
!     type(mineral_c), allocatable :: aux(:)
!     logical :: flag
!     aux=this%minerals
!     deallocate(this%minerals)
!     call this%set_num_mins_min_zone(this%num_minerals-size(old_min_ind))
!     call this%allocate_minerals_min_zone()
!     j=1
!     k=1
    
!     do 
!         i=1
!         flag=.true.
!         do
!             if (old_min_ind(i)==j) then
!                 flag=.false.
!                 exit
!             else if (i<size(old_min_ind)) then
!                 i=i+1
!             else
!                 exit
!             end if
!         end do
!         if (flag.eqv..true.) then
!             this%minerals(k)=aux(j)
!             if (k<this%num_minerals) then
!                 k=k+1
!             end if
!         end if
!         if (j<size(aux)) then
!             j=j+1
!         end if
!     end do
! end subroutine update_mineral_zone