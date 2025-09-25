subroutine compute_lumped_A_mat(this,A_mat_lumped)
!> A^L_i,i=A_i,i-1+A_i,i+A_i,i+1
!> A^L is the identity matrix if flow uniform
    use PDE_transient_m, only: PDE_1D_transient_c, diag_matrix_c
    implicit none
    
    class(PDE_1D_transient_c) :: this
    type(diag_matrix_c), intent(out) :: A_mat_lumped !> must be allocated
    
    integer(kind=4) :: i
    
    call A_mat_lumped%allocate_array(this%spatial_discr%Num_targets)
    
    A_mat_lumped%diag(1)=this%A_mat%diag(1)+this%A_mat%super(1)
    do i=2,this%spatial_discr%Num_targets-1
        A_mat_lumped%diag(i)=this%A_mat%sub(i-1)+this%A_mat%diag(i)+this%A_mat%super(i)
    end do
    A_mat_lumped%diag(this%spatial_discr%Num_targets)=this%A_mat%sub(this%spatial_discr%Num_targets-1)+&
    this%A_mat%diag(this%spatial_discr%Num_targets)
end subroutine