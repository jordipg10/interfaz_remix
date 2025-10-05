!> Computes mineral kinetic reaction rate
subroutine compute_rk_mineral(this,act_cat,saturation,react_surf,temp,rk)
    use kin_mineral_m, only: kin_mineral_c
    implicit none
    class(kin_mineral_c), intent(in) :: this
    real(kind=8), intent(in) :: act_cat(:) !> activities of catalysers (dim=num_cat)
    real(kind=8), intent(in) :: saturation !> saturation of reaction
    real(kind=8), intent(in) :: react_surf !> specific surface mineral
    real(kind=8), intent(in) :: temp !> temperature (Kelvin)
    real(kind=8), intent(out) :: rk !> kinetic reaction rate
    
    integer(kind=4) :: zeta,j,i
    real(kind=8) :: prod
    real(kind=8), parameter :: R=8.31446261815324 !> J/(mol*K)
    
    rk=0d0
    do j=1,this%params%num_par_reacts
        prod=1d0
        do i=1,this%params%num_cat
            prod=prod*act_cat(i)**this%params%p(j,i)
        end do
        rk=rk+prod*this%params%k(j)*abs(((saturation**this%params%theta(j))-1d0))**this%params%eta(j)
    end do
    if (saturation>1d0) then !> precipitation
        zeta=1
    else if (saturation<1d0) then !> dissolution
        zeta=-1
    end if
    rk=rk*react_surf*exp(-this%params%act_energy/(R*temp))
    rk=rk*zeta
end subroutine