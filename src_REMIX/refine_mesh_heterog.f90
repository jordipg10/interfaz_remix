subroutine refine_mesh_heterog(this,conc,conc_ext,rel_tol)
    use spatial_discr_1D_m, only: mesh_1D_Euler_heterog_c
    use vectors_m
    implicit none
    class(mesh_1D_Euler_heterog_c) :: this
    real(kind=8), intent(inout) :: conc(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(inout) :: conc_ext(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(in) :: rel_tol !> relative tolerance
    
    integer(kind=4) :: j,n
    real(kind=8), allocatable :: Delta_x_new(:)
    
    n=this%Num_targets
    
    do j=1,n-1
        if (inf_norm_vec_real(conc(:,j)-conc(:,j+1))/inf_norm_vec_real(conc(:,j))>=rel_tol) then
            Delta_x_new=[Delta_x_new,this%Delta_x(j)/2,this%Delta_x(j)/2]
        else
            Delta_x_new=[Delta_x_new,this%Delta_x(j)]
        end if
    end do
    deallocate(this%Delta_x)
    this%Delta_x=Delta_x_new
    this%Num_targets=size(this%Delta_x)
end subroutine