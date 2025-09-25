!subroutine initialise_target_waters(this,initial_water_types)
!    use chemistry_m
!    implicit none
!    class(chemistry_c) :: this
!    class(aqueous_chemistry_c), intent(in) :: initial_water_types(:)
!    
!    integer(kind=4) :: i,j,m,niter
!    real(kind=8), allocatable :: conc_init_species(:,:),u_init(:,:),c1_init(:,:),c2aq_init(:,:),gamma_1(:),gamma_2aq(:)
!    real(kind=8), parameter :: rel_tol=1d-3
!    integer(kind=4), parameter :: niter_max=30
!    
!    call this%allocate_target_waters(this%num_target_solids) !> suponemos biyeccion target waters - target solids
!!> Autentica chapuza ad hoc
!    if (size(initial_water_types)==2) then
!        do j=1,floor(this%num_target_waters/2d0)
!            this%target_waters_init(j)=initial_water_types(1)
!        end do
!        do j=1,ceiling(this%num_target_waters/2d0)
!            this%target_waters_init(j+floor(this%num_target_waters/2d0))=initial_water_types(2)
!        end do
!    end if
!    do j=1,this%num_target_waters
!        call this%target_waters_init(j)%set_solid_chemistry(this%target_solids(j))
!        !call this%target_waters_init(j)%set_conc_cst_act_species([1d0]) !> we assume we have ideal water
!        !call this%target_waters_init(j)%set_conc_sec_var_act_species() !> we initialise concentration secondary variable activity species
!        call this%target_waters_init(j)%compute_c2nc_aq_from_c1_aq_expl()
!        call this%target_waters_init(j)%compute_c2nc_from_c1_Picard(niter)
!        call this%target_waters_init(j)%compute_conc_comp_aq()
!    end do
!    this%target_waters=this%target_waters_init
!end subroutine