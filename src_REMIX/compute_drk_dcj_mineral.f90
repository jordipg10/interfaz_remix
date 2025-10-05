!subroutine compute_drk_dcj_mineral(this,j,c_j,rk_j,drk_dcj)
!    use kin_mineral_m
!    implicit none
!    class(kin_mineral_c), intent(in) :: this
!    integer(kind=4), intent(in) :: j
!    real(kind=8), intent(in) :: c_j
!    real(kind=8), intent(in) :: rk_j
!    real(kind=8), intent(out) :: drk_dcj
!        
!    drk_dcj=(rk_j*this%stoichiometry(j)*this%params%eta(j)*this%params%theta(j)*this%saturation**this%params%theta(j))/(c_j*(this%saturation**this%params%theta(j)-1d0))
!end subroutine