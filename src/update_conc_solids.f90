!> Updates concentration solids in Newton method
subroutine update_conc_solids(this,Delta_c_s,control_factor)
    use solid_chemistry_m
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_s(:) !> solid concentration difference
    real(kind=8), intent(in) :: control_factor !> must \f$\in (0,1)\f$
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (control_factor>1d0 .or. control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,size(this%concentrations)
        if (this%concentrations(i)+Delta_c_s(i)<=control_factor*this%concentrations(i)) then
            this%concentrations(i)=control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_s(i)>=this%concentrations(i)/control_factor) then
            this%concentrations(i)=this%concentrations(i)/control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_s(i)
        end if
        Delta_c_s(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine