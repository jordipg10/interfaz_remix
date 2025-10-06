subroutine refine_mesh_homog(this,conc,conc_ext,rel_tol)
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c
    use vectors_m
    implicit none
    class(mesh_1D_Euler_homog_c) :: this
    real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
    real(kind=8), intent(in) :: rel_tol !> relative tolerance
    
    real(kind=8), allocatable :: conc_ref(:,:)
    real(kind=8), allocatable :: conc_ext_ref(:,:)
    
    integer(kind=4) :: j,n,ratio
    real(kind=8) :: Delta_x_old,Delta_x_new
    
    n=this%Num_targets
    Delta_x_old=this%Delta_x
    Delta_x_new=Delta_x_old
    do j=1,n-1
        if (inf_norm_vec_real(conc(:,j)-conc(:,j+1))/inf_norm_vec_real(conc(:,j))>=rel_tol) then
            Delta_x_new=Delta_x_old/2
            Delta_x_old=Delta_x_new
        end if
    end do
    if (Delta_x_new<this%Delta_x) then
        ratio=nint(this%Delta_x/Delta_x_new)
        this%Num_targets=this%Num_targets*ratio
        this%Delta_x=Delta_x_new
        allocate(conc_ref(size(conc,1),this%Num_targets),conc_ext_ref(size(conc_ext,1),this%Num_targets))
        do j=1,this%Num_targets
            conc_ref(:,j)=conc(:,ceiling(j/(1d0*ratio)))
            conc_ext_ref(:,j)=conc_ext(:,ceiling(j/(1d0*ratio)))
        end do
        deallocate(conc,conc_ext)
        conc=conc_ref
        conc_ext=conc_ext_ref
    end if
end subroutine