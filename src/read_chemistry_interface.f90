!> \file read_chemistry_interface.f90
!> \brief Reads all chemistry input data from files.
!> \details
!> Dispatches to CHEPROO, PHREEQC, or PFLOTRAN readers depending
!> on the selected input format option.  Populates the chemistry
!> object with the chemical system, initial conditions, boundary
!> conditions, and reactive zone definitions.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup io
!> \see chemistry_m, read_chem_system_CHEPROO, read_chemistry_CHEPROO

!> This subroutine reads chemistry data depending on the selected option
subroutine read_chemistry_interface(this,root,path_pb,path_DB,num_tar)
    use chemistry_m, only: chemistry_c
    implicit none
    !> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root of input and output files of problem
    character(len=*), intent(in) :: path_pb !> path for input and output files of problem
    character(len=*), intent(in) :: path_DB !> path for chemical databases
    integer(kind=4), intent(in) :: num_tar !> expected number of target (domain) waters from the mesh

    !> Variables
    integer(kind=4) :: unit_chem_syst !> chemical system file unit
    integer(kind=4) :: unit_comp_opts !> computation options file unit
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    integer(kind=4) :: read_opt !> option for reading chemistry data
    
    !> File units (arbitrary)
    unit_chem_syst=1
    unit_comp_opts=51
    unit_loc_chem=2

    !> Read chemical data depending on the selected option (we impose CHEPROO here)
    read_opt=1 !> CHEPROO-based (chapuza)
    call this%set_read_opt(read_opt) !> set reading option (chapuza)
    if (this%read_opt==1) then !> CHEPROO-based
        call this%read_chemistry_CHEPROO(root,path_pb,path_DB,unit_chem_syst,unit_loc_chem,num_tar)
    else if (this%read_opt==2) then !> PHREEQC
        error stop "PHREEQC data input not fully implemented yet"
    else if (this%read_opt==3) then !> PFLOTRAN
        error stop "PFLOTRAN data input not fully implemented yet"
    else
        error stop "This data input option is not implemented yet"
    end if
end subroutine