!> Updates concentration aqueous species in Newton method
subroutine update_conc_aq_species(this,Delta_c_aq)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_aq(:)
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (this%solid_chemistry%reactive_zone%CV_params%control_factor>1d0 .or. &
    this%solid_chemistry%reactive_zone%CV_params%control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,this%aq_phase%num_species
        if (this%concentrations(i)+Delta_c_aq(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
        this%concentrations(i)) then
            this%concentrations(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_aq(i)>=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(i)=this%concentrations(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_aq(i)
        end if
        Delta_c_aq(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine update_conc_aq_species