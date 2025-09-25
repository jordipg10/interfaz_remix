!> Appends phase object to array of phase objects
subroutine append_phase(this,phase)
    use chem_system_m, only: chem_system_c, mineral_c, phase_c, gas_phase_c
    implicit none
    class(chem_system_c) :: this
    class(phase_c), intent(in) :: phase
    
    integer(kind=4) :: i
    type(mineral_c), allocatable :: mins(:)
    type(gas_phase_c), allocatable :: gas_phases(:)
    
    select type (phase)
    type is (mineral_c)
        mins=this%minerals
        if (allocated(this%minerals)) then
            deallocate(this%minerals)
        end if
        call this%allocate_minerals(size(mins)+1)
        do i=1,this%num_minerals-1
            this%minerals(i)=mins(i)
        end do
        this%minerals(this%num_minerals)=phase
    !type is (gas_phase_c)
    !>    gas_phases=this%gases
    !>    if (allocated(this%gas_phases)) then
    !>        deallocate(this%gas_phases)
    !>    end if
    !>    call this%allocate_gases(size(gases)+1)
    !>    do i=1,this%num_species-1
    !>        this%gas_phases(i)=gases(i)
    !>    end do
    !>    this%gases(this%num_species)=phase
    end select
end subroutine 