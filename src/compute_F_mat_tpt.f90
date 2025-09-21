subroutine compute_F_mat_tpt(this) !> diagonal matrix
!> F_ii=phi_i
    use transport_transient_m, only: transport_1D_transient_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    implicit none
    class(transport_1D_transient_c) :: this
    
    integer(kind=4) :: i,n
    real(kind=8) :: r_i !> radius (for radial symmetry)
    
    n=this%spatial_discr%Num_targets
     
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    
    
    this%F_mat%diag=this%tpt_props_heterog%porosity
    
    ! if (this%tpt_props_heterog%homog_flag .eqv. .true.) then
    !     this%F_mat%diag=this%tpt_props_heterog%porosity(1)
    ! else
    !select type (mesh=>this%spatial_discr)
    !type is (mesh_1D_Euler_homog_c)
    !    this%F_mat%diag=this%tpt_props_heterog%porosity(1)
    !type is (mesh_1D_Euler_heterog_c)
    !    this%F_mat%diag=this%tpt_props_heterog%porosity
    !type is (spatial_discr_rad_c)
    !    if (this%dimless.eqv..false.) then
    !        do i=1,n
    !            !r_i=mesh%targets(i)%coord(1) !> get radius
    !            this%F_mat%diag(i)=this%tpt_props_heterog%porosity(i)*mesh%Delta_r(i)*&
    !                (mesh%targets(i)%coord(1)**(mesh%dim-1))
    !        end do
    !    else
    !        do i=1,n
    !            !r_i=sum(mesh%Delta_r(1:i)) !> compute radius for radial symmetry
    !            this%F_mat%diag(i)=mesh%Delta_r(i)*(mesh%targets(i)%coord(1)**(mesh%dim-1))
    !        end do
    !    end if
    !end select
    !end if
    ! print *, this%tpt_props_heterog%porosity
    !print *, this%F_mat%diag
end subroutine