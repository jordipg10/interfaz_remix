!subroutine initialise_ext_waters(this)
!    use chemistry_m
!    implicit none
!    class(chemistry_c) :: this
!    
!    integer(kind=4) :: i,j
!    real(kind=8), allocatable :: c_ext(:)
!    
!    call this%set_num_ext_waters() !> suponemos biyeccion entre external waters y target solids 
!    call this%allocate_ext_waters()
!        
!    do j=1,this%num_ext_waters
!        call this%ext_waters(j)%set_solid_chemistry(this%target_solids(j))
!        call this%ext_waters(j)%set_aq_phase(this%aq_phase)
!        call this%ext_waters(j)%allocate_conc_aq_species()
!        this%ext_waters(j)%concentrations=0d0 !> no sink/source terms
!        !call this%ext_waters(j)%compute_conc_comp_aq()
!    end do
!end subroutine