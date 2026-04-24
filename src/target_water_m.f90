!> \file target_water_m.f90
!> \brief Target water class associating a spatial position with aqueous chemistry.
!> \details
!> Defines the target_water_c type used in Lagrangian and mixed
!> Eulerian-Lagrangian reactive transport formulations. Each target
!> water carries:
!>   - A unique identifier.
!>   - A spatial position vector (current and previous time level).
!>   - A pointer to its aqueous_chemistry_c object.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup transport
!> \see aqueous_chemistry_m, chemistry_m
module target_water_m
use aqueous_chemistry_m, only: aqueous_chemistry_c
implicit none
save
private
    type, public :: target_water_c !< Target water class
        integer(kind=4) :: id                           !< ID of this target water
        class(aqueous_chemistry_c), pointer :: aq_chem=>null()  !< Pointer to the aqueous chemistry class
        real(kind=8), allocatable :: pos(:)                       !< Spatial position of this water
        real(kind=8), allocatable :: pos_old(:)                       !< Spatial position of this water
    contains
        procedure :: set_id
        procedure :: set_pos
        procedure :: set_aq_chem
        procedure :: copy_tar_wat
        procedure :: update_old_attributes
    end type target_water_c

    contains
    
    subroutine set_pos(this, pos)
    class(target_water_c), intent(inout) :: this
    real(kind=8), intent(in) :: pos(:)
    this%pos = pos
    end subroutine set_pos
    
    subroutine set_aq_chem(this, aq_chem)
    class(target_water_c), intent(inout) :: this
    class(aqueous_chemistry_c), intent(in), target :: aq_chem
    if (.not. associated(aq_chem%aq_phase)) then
        error stop "Error in set_aq_chem: aq_phase must be associated."
    end if
    this%aq_chem => aq_chem
    end subroutine set_aq_chem
    
    subroutine copy_tar_wat(this, other)
    class(target_water_c), intent(inout) :: this
    class(target_water_c), intent(in), target :: other
    !type(aqueous_chemistry_c), target :: aq_chem_copy
    call this%set_id(other%id)
    call this%set_pos(other%pos)
    if (associated(other%aq_chem)) then
        !call aq_chem_copy%copy_aq_chem(other%aq_chem)
        !call this%aq_chem%copy_aq_chem(other%aq_chem)
        call this%set_aq_chem(other%aq_chem)
    else
        error stop "Error in copy_tar_wat: other%aq_chem is not associated."
    end if
    end subroutine copy_tar_wat
    
    subroutine set_id(this, id)
    class(target_water_c), intent(inout) :: this
    integer(kind=4), intent(in) :: id
    if (id.lt.0) then
        error stop "Error in set_id: id must be non-negative."
    end if
    this%id = id
    end subroutine set_id

    subroutine update_old_attributes(this)
    class(target_water_c), intent(inout) :: this
    ! if (allocated(this%pos)) then
    !     if (.not. allocated(this%pos_old)) then
    !         allocate(this%pos_old(size(this%pos)))
    !     end if
        this%pos_old = this%pos
    ! else
    !     error stop "Error in update_old_attributes: pos is not allocated."
    ! end if
    call this%aq_chem%update_old_attributes()
    end subroutine update_old_attributes

end module target_water_m
