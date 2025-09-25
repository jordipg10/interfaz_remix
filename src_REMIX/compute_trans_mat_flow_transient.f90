subroutine compute_trans_mat_flow_transient(this)
    !use PDE_m, only: PDE_1D_c
    !use transport_m, only: transport_1D_c, diffusion_1D_c
    use flow_transient_m, only: flow_transient_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c, spatial_discr_c
    implicit none
    class(flow_transient_c) :: this
    real(kind=8) :: r_i,h_i,r_min_D
    integer(kind=4) :: i,n,opcion
    
    n=this%spatial_discr%Num_targets

    if (this%spatial_discr%adapt_ref.eq.1) then
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super)
        call this%allocate_trans_mat()
    end if
    select type (mesh=>this%spatial_discr)
    type is (spatial_discr_rad_c)
        if (this%dimless.eqv..true.) then
            !r_min_D=mesh%r_min/this%char_params_flow%char_length
            do i=1,n-1
                r_i=mesh%r_min_D+sum(mesh%Delta_r_D(1:i)) !> $r_D,i$ radial coordinate
                h_i=(mesh%targets(i+1)%coord_D(1)-mesh%targets(i)%coord_D(1)) !> $h_D,i$ radial step
                this%trans_mat%sub(i)=2d0*(r_i**(mesh%dim-1))/h_i !> sub-diagonal element
            end do
        end if
    end select
    this%trans_mat%super=this%trans_mat%sub
    !this%trans_mat%diag=0d0 !> initialise diagonal elements
    this%trans_mat%diag(2:n-1)=-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1)
end subroutine