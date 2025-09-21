!subroutine compute_A_mat_ODE_tpt(this,A_mat)
!!> dy/dt=-Ay+b
!!> A=-F^(-1)*T: tridiagonal, symmetric, positive semi-definite
!>    use transport_transient_1_m
!>    !use matrices_m
!>    implicit none
!>    class(transport_1D_transient_c), intent(in) :: this
!>    class(array_c), pointer, intent(out) :: A_mat
!>    !real(kind=8), intent(out), optional :: A_sub(:),A_diag(:),A_super(:)
!>    
!>    integer(kind=4) :: n
!>    !real(kind=8), allocatable :: F_mat(:),G_mat(:,:),inv(:,:),aux_mat(:,:)
!>    type(tridiag_matrix_c), target :: A_mat_tpt
!>    !select type (this)
!>    !class is (diffusion_transient_c)
!>        n=this%spatial_discr%Num_targets
!>        !allocate(G_mat(n,n),inv(n,n),aux_mat(n,n))
!>        !A_mat=0d0
!>        !A_mat(1,1:2)=-this%trans_mat(1,1:2)
!>        !forall (i=2:n-1)
!>        !>    A_mat(i,i-1)=-this%trans_mat(i,i-1)
!>        !>    A_mat(i,i)=-this%trans_mat(i,i)
!>        !>    A_mat(i,i+1)=-this%trans_mat(i,i+1)
!>        !end forall
!>        !A_mat(n,n-1:n)=-this%trans_mat(n,n-1:n)
!>        !F_mat=this%compute_F_mat()
!>        !G_mat=0d0
!>        !G_mat(1,1:2)=F_mat(1)
!>        !A_mat(1,:)=A_mat(1,:)/F_mat(1)
!>        !A_mat(1,1)=A_mat(1,1)/F_mat(1)
!>        !forall (i=2:n-1)
!>        !>    G_mat(i,i-1)=F_mat(i)
!>        !>    G_mat(i,i)=F_mat(i)
!>        !>    G_mat(i,i+1)=F_mat(i)
!>        !>    !A_mat(i,:)=A_mat(i,:)/F_mat(i)
!>        !>    A_mat(i,i)=A_mat(i,i)/F_mat(i)
!>        !end forall
!>        !A_mat(n,:)=A_mat(n,:)/F_mat(n)
!>        !A_mat(n,n)=A_mat(n,n)/F_mat(n)
!>        !G_mat(n,n-1:n)=F_mat(n)
!>        !call inv_matrix(sqrt(G_mat),inv)
!>        !print *, det(inv)
!>        !aux_mat=matmul(A_mat,inv)
!>        !A_mat=matmul(sqrt(G_mat),aux_mat)
!>        
!>        select type (trans_mat=>this%trans_mat)
!>        type is (tridiag_matrix_c)
!>            A_mat_tpt%diag=trans_mat%diag/this%F_mat%diag
!>            A_mat_tpt%sub=trans_mat%sub
!>            A_mat_tpt%super=trans_mat%super
!>            A_mat=>A_mat_tpt
!>        end select
!>        
!>        !if (present(A_sub) .and. present(A_diag) .and. present(A_super)) then
!>        !>    forall (i=1:n-1)
!>        !>        A_sub(i)=A_mat(i+1,i)
!>        !>        A_super(i)=A_mat(i,i+1)
!>        !>        A_diag(i)=A_mat(i,i)
!>        !>    end forall
!>        !>    A_diag(n)=A_mat(n,n)
!>        !end if
!>    !end select
!end subroutine