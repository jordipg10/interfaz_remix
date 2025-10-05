!!> Gets the indices of the subset of target waters that will mix with the target water number tar_wat_indices(i)
!subroutine get_mixing_waters_indices(i,tar_wat_indices,mix_wat_indices)
!    use chemistry_m
!    implicit none
!    integer(kind=4), intent(in) :: i
!    integer(kind=4), intent(in) :: tar_wat_indices(:)
!    integer(kind=4), intent(out) :: mix_wat_indices(:) ! already allocated
!    
!    integer(kind=4) :: j,k
!    
!    !if (i==1) then
!    !    mix_wat_indices=[tar_wat_indices(2)]
!    !else if (i==size(tar_wat_indices)) then
!    !    mix_wat_indices=[tar_wat_indices(size(tar_wat_indices)-1)]
!    !else
!    !    mix_wat_indices=[tar_wat_indices(i-1),tar_wat_indices(i+1)]
!    !end if
!    !k=0
!    !do j=1,size(tar_wat_indices)
!    !    if (j/=i) then
!    !        k=k+1
!    !        mix_wat_indices(k)=tar_wat_indices(j)
!    !    end if
!    !end do
!    !mix_wat_indices=this%mixing_waters_indices%cols(tar_wat_indices(i))%col_1
!end subroutine
!        