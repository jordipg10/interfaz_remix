!> Updates concentration aqueous primary species in iterative method
subroutine update_conc_aq_prim_species(this,Delta_c1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq
    real(kind=8), allocatable :: c1_old(:)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    
    if (n_p_aq/=size(Delta_c1) .and. this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species/=size(Delta_c1)) then
        error stop "Dimension error in update_conc_aq_prim_species"
    end if
    c1_old=this%concentrations(this%indices_aq_phase(1:n_p_aq))
    do i=1,n_p_aq
        if (this%concentrations(this%indices_aq_phase(i))+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
        this%concentrations(this%indices_aq_phase(i))) then
            this%concentrations(this%indices_aq_phase(i))=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(this%indices_aq_phase(i))
        else if (this%concentrations(this%indices_aq_phase(i))+Delta_c1(i)>=this%concentrations(this%indices_aq_phase(i))/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(this%indices_aq_phase(i))=this%concentrations(this%indices_aq_phase(i))/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(this%indices_aq_phase(i))=this%concentrations(this%indices_aq_phase(i))+Delta_c1(i)
        end if
        Delta_c1(i)=this%concentrations(this%indices_aq_phase(i))-c1_old(i)
    end do
end subroutine