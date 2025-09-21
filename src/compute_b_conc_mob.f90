!> Computes vector b in Euler implicit linear system for mobile zone concentrations in MRMT
!> A*c_mob^(k+1)=b
subroutine compute_b_conc_mob(this,theta,Delta_t,conc_mob_old,conc_imm_old,b)
    use MRMT_m, only: MRMT_c
    use matrices_m, only: prod_mat_vec,tridiag_matrix_c
    implicit none
    class(MRMT_c), intent(in) :: this
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: conc_mob_old(:) !> c_mob^k
    real(kind=8), intent(in) :: conc_imm_old(:) !> c_imm^k
    real(kind=8), intent(out) :: b(:) 
    
    integer(kind=8) :: i,j,n_imm
    real(kind=8) :: sum1
    type(tridiag_matrix_c) :: X_mat
    real(kind=8), allocatable :: f(:),sum2(:)
    
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1"
    X_mat%sub=((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%sub
    X_mat%super=((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%super
    sum1=0d0
    do i=1,this%n_imm
        sum1=sum1+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(1d0-theta-(&
        theta*this%imm_zones(i)%exch_rate*Delta_t*(1d0-theta))/(1d0+this%imm_zones(i)%exch_rate*Delta_t*theta))
    end do
    X_mat%diag=1d0+((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%diag-(Delta_t/this%mob_zone%mob_por)*sum1
    allocate(sum2(size(conc_imm_old)))
    sum2=0d0
    do i=1,this%n_imm
        sum2=sum2+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(theta-1-theta*(&
        1-this%imm_zones(i)%exch_rate*Delta_t*(1-theta))/(1+this%imm_zones(i)%exch_rate*Delta_t*theta))*conc_imm_old(i)
    end do
    f=(Delta_t/this%mob_zone%mob_por)*sum2
    b=prod_mat_vec(X_mat,conc_mob_old)-f
end subroutine