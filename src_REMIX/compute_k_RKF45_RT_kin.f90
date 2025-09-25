!!> f(c)=L(c)*inv(F)+transp(S_k)*rk(c)*inv(F)=(c*transp(T)+c_ext+transp(S_k)*rk(c))*inv(F)
!function compute_k_RKF45_RT_kin(this,F_mat,trans_mat,Delta_t,conc_RK4,conc_ext) result(k)
!>    use chemistry_kin_m
!>    !use vectores_m
!>    implicit none
!>    class(chemistry_kin_c), intent(in) :: this
!>    !integer(kind=4), intent(in) :: j !> target
!>    class(diag_matrix_c), intent(in) :: F_mat
!>    class(tridiag_matrix_c), intent(in) :: trans_mat
!>    !real(kind=8), intent(in) :: source_term_PDE(:,:)
!>    real(kind=8), intent(in) :: Delta_t
!>    real(kind=8), intent(in) :: conc_RK4(:,:)
!>    real(kind=8), intent(in) :: conc_ext(:,:)
!>    real(kind=8), allocatable :: k(:,:,:)
!>    
!>    integer(kind=4) :: n,i,icol,Num_output
!>    real(kind=8) :: Delta_t_old,Delta_t_new,Time
!>    real(kind=8), allocatable :: conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:),aux_mat(:,:),x_mat(:,:)
!>    
!>    !allocate(k(size(conc_RK4,1),6))
!>    allocate(k(size(conc_RK4,1),size(conc_RK4,2),6))
!>    !allocate(x_mat(size(conc_RK4,1),size(conc_RK4,2)))
!>    
!>    !k(:,1)=this%compute_f_RKF45_RT_kin(j,conc_RK4,conc_ext,F_mat,trans_mat)
!>    k(:,:,1)=this%compute_f_RKF45_RT_kin(conc_RK4,conc_ext,F_mat,trans_mat)
!>    !k(:,2)=this%compute_f_RKF45_RT_kin(j,conc_RK4+Delta_t*k(:,1)/4,conc_ext,F_mat,trans_mat)
!>    k(:,:,2)=this%compute_f_RKF45_RT_kin(conc_RK4+Delta_t*k(:,:,1)/4,conc_ext,F_mat,trans_mat)
!>    !!k(:,j,2)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*k(:,1)/4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !k(:,3)=this%compute_f_RKF45_RT_kin(j,conc_RK4+Delta_t*(3*k(:,1)+9*k(:,2))/32,conc_ext,F_mat,trans_mat)
!>    k(:,:,3)=this%compute_f_RKF45_RT_kin(conc_RK4+Delta_t*(3*k(:,:,1)+9*k(:,:,2))/32,conc_ext,F_mat,trans_mat)
!>    !!k(:,j,3)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(3*k(:,1)+9*k(:,2))/32)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !k(:,4)=this%compute_f_RKF45_RT_kin(j,conc_RK4+Delta_t*(1932*k(:,1)-7200*k(:,2)+7296*k(:,3))/2197,conc_ext,F_mat,trans_mat)
!>    k(:,:,4)=this%compute_f_RKF45_RT_kin(conc_RK4+Delta_t*(1932*k(:,:,1)-7200*k(:,:,2)+7296*k(:,:,3))/2197,conc_ext,F_mat,trans_mat)    
!>    !!k(:,j,4)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(1932*k(:,1)-7200*k(:,2)+7296*k(:,3))/2197)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !k(:,5)=this%compute_f_RKF45_RT_kin(j,conc_RK4+Delta_t*(439*k(:,1)/216-8*k(:,2)+3680*k(:,3)/513-845*k(:,4))/4104,conc_ext,F_mat,trans_mat)
!>    k(:,:,5)=this%compute_f_RKF45_RT_kin(conc_RK4+Delta_t*(439*k(:,:,1)/216-8*k(:,:,2)+3680*k(:,:,3)/513-845*k(:,:,4))/4104,conc_ext,F_mat,trans_mat)    
!>    !!k(:,j,5)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(439*k(:,1)/216-8*k(:,2)+3680*k(:,3)/513-845*k(:,4)/4104))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !k(:,6)=this%compute_f_RKF45_RT_kin(j,conc_RK4+Delta_t*(-8*k(:,1)/27+2*k(:,2)-3544*k(:,3)/2565+1859*k(:,4)/4104-11*k(:,5))/40,conc_ext,F_mat,trans_mat)
!>    k(:,:,6)=this%compute_f_RKF45_RT_kin(conc_RK4+Delta_t*(-8*k(:,:,1)/27+2*k(:,:,2)-3544*k(:,:,3)/2565+1859*k(:,:,4)/4104-11*k(:,:,5))/40,conc_ext,F_mat,trans_mat)    
!>    !!k(:,j,6)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(-8*k(:,1)/27+2*k(:,2)-3544*k(:,3)/2565+1859*k(:,4)/4104-11*k(:,5)/40))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!
!>    !do j=1,n
!>    !>    k(:,j,1)=compute_f_RKF45_RT_kin(this,conc,conc_ext,F_mat,trans_mat
!>    !>    k(:,j,2)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*k(:,1)/4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !>    k(:,j,3)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(3*k(:,1)+9*k(:,2))/32)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !>    k(:,j,4)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(1932*k(:,1)-7200*k(:,2)+7296*k(:,3))/2197)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !>    k(:,j,5)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(439*k(:,1)/216-8*k(:,2)+3680*k(:,3)/513-845*k(:,4)/4104))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !>    k(:,j,6)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(-8*k(:,1)/27+2*k(:,2)-3544*k(:,3)/2565+1859*k(:,4)/4104-11*k(:,5)/40))/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
!>    !end do
!end function