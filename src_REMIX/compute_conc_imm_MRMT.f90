!> Computes concentration immobile zones at time k+1
subroutine compute_conc_imm_MRMT(this,theta,conc_imm_old,conc_mob_old,conc_mob_new,Delta_t,conc_imm_new)
    use MRMT_m, only: MRMT_c
    implicit none
    class(MRMT_c), intent(in) :: this
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: conc_imm_old(:) !> c_imm^k
    real(kind=8), intent(in) :: conc_mob_old(:) !> c_m^k
    real(kind=8), intent(in) :: conc_mob_new(:) !> c_m^(k+1)
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_imm_new(:) !> c_imm^(k+1)
    
    integer(kind=8) :: j
    
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1"
    do j=1,this%n_imm
        conc_imm_new(j)=(conc_imm_old(j)*(1d0-this%imm_zones(j)%exch_rate*Delta_t*(1-theta))+this%imm_zones(j)%exch_rate*Delta_t*&
        (theta*conc_mob_new(j)+(1-theta)*conc_mob_old(j)))/(1+this%imm_zones(j)%exch_rate*Delta_t*theta)    
    end do
end subroutine