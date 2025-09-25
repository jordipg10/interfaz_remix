!> Computes matrix A in Euler implicit linear system for mobile zone concentrations in MRMT
!> A*c_mob^(k+1)=b
subroutine compute_A_mat_conc_mob(this,theta,Delta_t,A_mat)
    use MRMT_m, only: MRMT_c
    use matrices_m, only: tridiag_matrix_c
    implicit none
    class(MRMT_c), intent(in) :: this
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    class(tridiag_matrix_c), intent(out) :: A_mat
    
    integer(kind=8) :: i,j
    real(kind=8) :: sum
    
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1"
    A_mat%sub=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%sub
    A_mat%super=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%super
    sum=0d0
    do i=1,this%n_imm
        sum=sum+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(theta-(this%imm_zones(i)%exch_rate*&
        Delta_t*theta**2)/(1+this%imm_zones(i)%exch_rate*Delta_t*theta))
    end do
    A_mat%diag=1d0-(theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%diag+(Delta_t/this%mob_zone%mob_por)*sum
end subroutine