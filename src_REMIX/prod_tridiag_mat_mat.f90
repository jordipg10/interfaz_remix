function prod_tridiag_mat_mat(this,B_mat) result(C_mat) !> AB=C
    use matrices_m, only : tridiag_matrix_c
    implicit none
    class(tridiag_matrix_c), intent(in) :: this !> nxn
    real(kind=8), intent(in) :: B_mat(:,:) !> nxm
    real(kind=8), allocatable :: C_mat(:,:) !> nxm
    
    integer(kind=4) :: i,j,n,m
    
    m=size(B_mat,2)
    !select type (this)
    !type is (tridiag_matrix_c)
        n=size(this%diag)
        if (n/=size(B_mat,1)) error stop "Dimension error in prod_mat_mat"
        if (.not. allocated(C_mat)) then
            allocate(C_mat(n,m))
        end if
        C_mat(1,:)=this%diag(1)*B_mat(1,:)+this%super(1)*B_mat(2,:)
        do i=2,n-1
            do j=1,m
                C_mat(i,j)=this%sub(i-1)*B_mat(i-1,j)+this%diag(i)*B_mat(i,j)+this%super(i)*B_mat(i+1,j)
            end do
        end do
        C_mat(n,:)=this%sub(n-1)*B_mat(n-1,:)+this%diag(n)*B_mat(n,:)
    !end select
end function 