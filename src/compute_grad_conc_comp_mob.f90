!subroutine compute_grad_conc_comp(this,react_zone_ind,Delta_x,BCs,grad_conc_comp)
!    use chemistry_m
!    use BCs_m
!    implicit none
!    class(chemistry_c), intent(in) :: this
!    integer(kind=4), intent(in) :: react_zone_ind ! reactive zone index
!    real(kind=8), intent(in) :: Delta_x
!    class(BCs_t), intent(in) :: BCs
!    real(kind=8), intent(out) :: grad_conc_comp(:,:) ! must be already allocated
!    
!    integer(kind=4) :: i,j,num_tar_wat
!    integer(kind=4), allocatable :: tw_indices(:)
!    
!    call this%link_target_waters_reactive_zone(react_zone_ind,tw_indices)
!    num_tar_wat=size(tw_indices)
!    
!    if (size(grad_conc_comp,1)/=this%target_waters(tw_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species .or. size(grad_conc_comp,2)/=num_tar_wat) error stop " Dimension error in compute_grad_conc_comp"
!    
!    if (BCs%BCs_label(1)==2 .and. BCs%BCs_label(2)==2) then ! Neumann homogeneous
!        do i=1,this%target_waters(tw_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
!            grad_conc_comp(i,1)=0d0
!            do j=2,num_tar_wat-1
!                grad_conc_comp(i,j)=(this%target_waters(tw_indices(j+1))%conc_comp(i)-this%target_waters(tw_indices(j-1))%conc_comp(i))/(2d0*Delta_x) ! CFD
!            end do
!            grad_conc_comp(i,num_tar_wat)=0d0
!        end do
!    end if
!end subroutine