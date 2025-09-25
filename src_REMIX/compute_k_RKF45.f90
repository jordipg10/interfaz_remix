 function compute_k_RKF45(this,Delta_t,conc_RK4) result(k)
    use PDE_transient_m, only: PDE_1D_transient_c, prod_mat_vec
    implicit none
    class(PDE_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: Delta_t
    real(kind=8), intent(in) :: conc_RK4(:)
    real(kind=8), allocatable :: k(:,:)
    
    integer(kind=4) :: n,i,j,icol,Num_output
    real(kind=8) :: Delta_t_old,Delta_t_new,Time
    real(kind=8), allocatable :: conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:)
    
    allocate(k(size(conc_RK4),6))
    
    k(:,1)=prod_mat_vec(this%trans_mat,conc_RK4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
    k(:,2)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*k(:,1)/4)/this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
    k(:,3)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(3*k(:,1)+9*k(:,2))/32)/this%F_mat%diag +&
     this%source_term_PDE/this%F_mat%diag
    k(:,4)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(1932*k(:,1)-7200*k(:,2)+7296*k(:,3))/2197)/this%F_mat%diag + &
    this%source_term_PDE/this%F_mat%diag
    k(:,5)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(439*k(:,1)/216-8*k(:,2)+3680*k(:,3)/513-845*k(:,4)/4104))/this%F_mat%diag+&
     this%source_term_PDE/this%F_mat%diag
    k(:,6)=prod_mat_vec(this%trans_mat,conc_RK4+Delta_t*(-8*k(:,1)/27+2*k(:,2)-3544*k(:,3)/2565+1859*k(:,4)/4104-11*k(:,5)/40))/&
    this%F_mat%diag + this%source_term_PDE/this%F_mat%diag
    
end function