subroutine compute_F_mat_flow(this)
    use flow_transient_m, only: flow_transient_c, spatial_discr_rad_c
    implicit none
    class(flow_transient_c) :: this
    
    integer(kind=4) :: i,n
    real(kind=8) :: r_i, h_i, r_1, r_n, h_1, h_n

    !> Compute the F matrix for the flow transient object
    !> This subroutine computes the F matrix for the flow transient object based on its properties

    n=this%spatial_discr%Num_targets
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    select type (mesh=>this%spatial_discr)
    type is (spatial_discr_rad_c)
        if (this%dimless) then
            r_1=mesh%targets(1)%coord_D(1) !> $r_D,1$ radial coordinate
            h_1=mesh%targets(2)%coord_D(1)-mesh%r_min_D !> $h_D,1$ radial step
            this%F_mat%diag(1)=h_1*(r_1**(mesh%dim-1))
            do i=2,n-1
                r_i=mesh%targets(i)%coord_D(1) !> $r_D,i$ radial coordinate
                h_i=mesh%targets(i+1)%coord_D(1)-mesh%targets(i-1)%coord_D(1) !> $h_{D,i+1}-h_{D,i-1}$
                !r_plus=mesh%r_min+sum(mesh%Delta_r(1:i))-5d-1*mesh%Delta_r(i)
                !r_minus=mesh%r_min+sum(mesh%Delta_r(1:i))-5d-1*mesh%Delta_r(i)
                this%F_mat%diag(i)=h_i*(r_i**(mesh%dim-1))
            end do
            r_n=mesh%targets(n)%coord_D(1) !> $r_D,n$ radial coordinate
            h_n=mesh%r_max_D-mesh%targets(n-1)%coord_D(1) !> $h_D,n$ radial step
            this%F_mat%diag(n)=h_n*(r_n**(mesh%dim-1))
        else
            do i=1,n
                r_i=mesh%r_min+sum(mesh%Delta_r(1:i))-5d-1*mesh%Delta_r(i)
                this%F_mat%diag(i)=mesh%Delta_r(i)*(r_i**(mesh%dim-1))*this%flow_props_heterog%storativity(i)
            end do
        end if
    end select
    !print *, "F_mat", this%F_mat%diag
end subroutine compute_F_mat_flow