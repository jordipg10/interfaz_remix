subroutine refine_mesh_rad(this,conc,conc_ext,rel_tol)
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    implicit none
    class(spatial_discr_rad_c) :: this
    real(kind=8), intent(inout) :: conc(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(inout) :: conc_ext(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(in) :: rel_tol !> relative tolerance
end subroutine