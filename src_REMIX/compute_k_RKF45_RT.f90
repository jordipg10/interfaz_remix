!!> x=U_m*c_nc*F
!!> f(x)=U_m*L(c_nc)+U_m*transp(S_k,nc)*rk(c_nc)=U_m*c_nc*transp(T)+U_m*g+U_m*transp(S_k,nc)*rk(c_nc)=x*inv(F)*transp(T)+U_m*g+U_m*transp(S_k,nc)*rk(c_nc)
!function compute_k_RKF45_RT(this,trans_mat,Delta_t,conc_RK4) result(k)
!>    use chemistry_eq_kin_m
!>    !use vectores_m
!>    implicit none
!>    class(chemistry_eq_kin_c), intent(in) :: this
!>    class(tridiag_matrix_c), intent(in) :: trans_mat
!>    !real(kind=8), intent(in) :: source_term_PDE(:,:)
!>    real(kind=8), intent(in) :: Delta_t
!>    real(kind=8), intent(in) :: conc_RK4(:,:)
!>    real(kind=8), allocatable :: k(:,:,:)
!>    
!>    integer(kind=4) :: n,i,j,icol,Num_output
!>    real(kind=8) :: Delta_t_old,Delta_t_new,Time
!>    real(kind=8), allocatable :: conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:),aux_mat(:,:)
!>    
!>    allocate(k(size(conc_RK4,1),size(conc_RK4,2),6))
!>    allocate(aux_mat(size(conc_RK4,1),size(conc_RK4,2)))
!>    
!>    
!>    !select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>    !type is (chem_system_eq_kin_c)
!>        do i=1,size(conc_RK4,1)
!>            aux_mat(i,:)=prod_mat_vec(trans_mat,conc_RK4(i,:))
!>            !aux_vec(i,1)=prod_mat_vec(trans_mat,conc_RK4(i,:))
!>        end do
!>        k(:,:,1)=matmul(this%solid_chemistry%reactive_zone%chem_syst%mob_comp_mat,aux_mat)
!>        !do j=1,n
!>        !>    k(:,j,1)=prod_mat_vec(chem_syst%mob_comp_mat,conc_RK4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !>    k(:,j,2)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*k(:,1)/4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !>    k(:,j,3)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(3*k(:,1)+9*k(:,2))/32)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !>    k(:,j,4)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(1932*k(:,1)-7200*k(:,2)+7296*k(:,3))/2197)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !>    k(:,j,5)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(439*k(:,1)/216-8*k(:,2)+3680*k(:,3)/513-845*k(:,4)/4104))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !>    k(:,j,6)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(-8*k(:,1)/27+2*k(:,2)-3544*k(:,3)/2565+1859*k(:,4)/4104-11*k(:,5)/40))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>        !end do
!end function