subroutine initialise_target_solids(this,n)
    use chemistry_m, only : chemistry_c
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: n !> number of target solids ( we assume = number targets transport )
    
    integer(kind=4) :: i,j
    real(kind=8), allocatable :: conc_init_species(:,:),u_init(:,:),c1_init(:,:),c2aq_init(:,:),gamma_1(:),gamma_2aq(:)
    
    call this%allocate_target_solids(n)
    
    do j=1,this%num_target_solids
        call this%target_solids(j)%set_reactive_zone(this%reactive_zones(1)) !> we assume only 1 reactive zone
        call this%target_solids(j)%allocate_conc_solids()
    end do
!> Totalmente arbitrario
    !do j=1,floor(this%num_target_solids/2d0)
    !>    this%target_solids(j)%concentrations=1d1
    !end do
    !do j=1,ceiling(this%num_target_solids/2d0)
    !>    this%target_solids(j+floor(this%num_target_solids/2d0))%concentrations=5d0
    !end do
    do j=1,this%num_target_solids
        this%target_solids(j)%concentrations=1d0
    end do
    !print *, this%reactive_zones(1)%non_flowing_species(1)%name
    
    !print *, associated(this%target_solids(1)%reactive_zone)
    !print *, this%target_solids(1)%reactive_zone%num_non_flowing_species
    !print *, this%target_solids(1)%reactive_zone%non_flowing_species(1)%name
end subroutine