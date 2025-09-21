subroutine compute_F_mat_diff(this) !> diagonal matrix
!> F_ii=phi_i*(i-1/2)^(d-1)
    use diffusion_transient_m, only: diffusion_1D_transient_c, spatial_discr_rad_c, diag_matrix_c
    implicit none
    class(diffusion_1D_transient_c) :: this
    
    integer(kind=4) :: i,n
    real(kind=8) :: r_i
    
    n=this%spatial_discr%Num_targets
    if (this%spatial_discr%adapt_ref==1) then
        deallocate(this%F_mat%diag)
        call this%F_mat%allocate_array(n)
    end if
    select type (mesh=>this%spatial_discr)
    type is (spatial_discr_rad_c)
        if (mesh%dim==1) then
            this%F_mat%diag=this%diff_props_heterog%porosity
        else
            forall (i=1:n)
                this%F_mat%diag(i)=this%diff_props_heterog%porosity(i)*(i-5d-1)**(mesh%dim-1)
            end forall
        end if
    end select
end subroutine