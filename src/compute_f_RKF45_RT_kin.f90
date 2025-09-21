!!> f(c)=L(c)*inv(F)+transp(S_k)*rk(c)*inv(F)=(c*transp(T)+r*c_ext+transp(S_k)*rk(c))*inv(F)
!function compute_f_RKF45_RT_kin(this,conc,conc_ext,F_mat,trans_mat) result(f)
!>    use chemistry_kin_m
!>    !use vectores_m
!>    implicit none
!>    class(chemistry_kin_c), intent(in) :: this
!>    !integer(kind=4), intent(in) :: j !> target
!>    real(kind=8), intent(in) :: conc(:,:)
!>    real(kind=8), intent(in) :: conc_ext(:,:)
!>    class(diag_matrix_c), intent(in) :: F_mat
!>    class(tridiag_matrix_c), intent(in) :: trans_mat !> T
!>    !real(kind=8), intent(in) :: source_term_PDE(:) !> g
!>    !real(kind=8), intent(in) :: conc(:,:)
!>    real(kind=8), allocatable :: f(:,:) !> tiene que estar alocatado
!>    
!>    integer(kind=4) :: n,i,icol,Num_output
!>    real(kind=8) :: Delta_t_old,Delta_t_new,Time
!>    real(kind=8), allocatable :: conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:),rk_vec(:),rk_mat(:,:)
!>    
!>    !allocate(k(size(conc_RK4,1),size(conc_RK4,2),6))
!>    !allocate(x_mat(size(conc_RK4,1),size(conc_RK4,2)))
!>    
!>    f=trans_mat%prod_mat_mat(transpose(conc))
!>    f=transpose(f)
!>    
!>    allocate(rk_mat(this%solid_chemistry%reactive_zone%chem_syst%num_reactions,size(F_mat%diag)))
!>    call this%compute_rk_mat(this%solid_chemistry%reactive_zone%chem_syst%species,this%solid_chemistry%reactive_zone%chem_syst%kin_reactions,conc,rk_mat)
!>    f=f+conc_ext+matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%stoich_mat),rk_mat)
!>    
!>    !do j=1,n
!>    !>    f(:,j)=f(:,j)/F_mat%diag(j)
!>    !end do
!>    
!>    !n=size(F_mat%diag)
!>    !
!>    !if (j.eqv.1) then
!>    !>    !do i=1,size(conc,1)
!>    !>        f=conc(:,j)*trans_mat%diag(j)+conc(:,j+1)*trans_mat%super(j)
!>    !else if (j.eqv.n) then
!>    !>    f=conc(:,j-1)*trans_mat%sub(j-1)+conc(:,j)*trans_mat%diag(j)
!>    !else
!>    !>    f=conc(:,j-1)*trans_mat%sub(j-1)+conc(:,j)*trans_mat%diag(j)+conc(:,j+1)*trans_mat%super(j)
!>    !end if
!>    !
!>    !allocate(rk_vec(this%solid_chemistry%reactive_zone%chem_syst%num_reactions))
!>    !call this%compute_rk_vec(j,this%solid_chemistry%reactive_zone%chem_syst%species,this%solid_chemistry%reactive_zone%chem_syst%kin_reactions,conc(:,j),rk_vec)
!>    !f=f+conc_ext+matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%stoich_mat),rk_vec)
!end function