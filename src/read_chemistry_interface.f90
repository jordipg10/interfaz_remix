!> \file read_chemistry_interface.f90
!> \brief Reads all chemistry input data for the reactive-mixing interface executable.
!> \details
!> Dispatches to CHEPROO, PHREEQC, or PFLOTRAN readers depending
!> on the selected input format option.  Populates the chemistry
!> object with the chemical system, initial conditions, boundary
!> conditions, and reactive zone definitions, and links the target
!> waters to their reactive zones.
!>
!> Unlike read_chemistry (used by the full REMIX program), this reader
!>   - does NOT read the computation options ("<root>_comp_opts.dat") nor the
!>     output options ("<root>_out_opts.dat"): the interface routines write
!>     their own output files and do not use those options.
!> The convergence parameters file "<root>_CV_params.dat" is optional (see
!> read_CV_params): without it the defaults of CV_params_s apply.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup io
!> \see chemistry_m, read_chemistry, read_chem_system_CHEPROO, read_chemistry_CHEPROO

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
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    integer(kind=4) :: read_opt !> option for reading chemistry data
    integer(kind=4) :: i !> reactive zone loop index
    
    !> File units (arbitrary)
    unit_chem_syst=1
    unit_loc_chem=2

    !> Convergence parameters ("<root>_CV_params.dat" is optional: defaults of CV_params_s otherwise).
    !> Read BEFORE the chemistry, since the reactive zones are linked to this%CV_params while they are read.
    call this%CV_params%read_CV_params(path_pb,root)
    !> Read chemical data depending on the selected option (we impose CHEPROO here)
    read_opt=1 !> CHEPROO-based (chapuza)
    call this%set_read_opt(read_opt) !> set reading option (chapuza)
    if (this%read_opt==1) then !> CHEPROO-based
        call this%read_chemistry_CHEPROO(root,path_pb,path_DB,unit_chem_syst,unit_loc_chem)
    else if (this%read_opt==2) then !> PHREEQC
        error stop "PHREEQC data input not fully implemented yet"
    else if (this%read_opt==3) then !> PFLOTRAN
        error stop "PFLOTRAN data input not fully implemented yet"
    else
        error stop "This data input option is not implemented yet"
    end if
    !> Check the number of domain (target) waters against the mesh
    if (this%num_target_waters/=num_tar) then
        write(*,'(A,I0,A,I0,A)') "WARNING: number of target waters read (", this%num_target_waters, &
            ") differs from the number of targets given (", num_tar, ")"
    end if
    !> Assign each reactive zone its canonical id (position in reactive_zones) BEFORE linking,
    !> so set_wat_indices_zones is called with the reactive-zone ids already set.
    do i=1,this%num_reactive_zones
        this%reactive_zones(i)%id=i
    end do
    call this%set_wat_indices_zones() !> link waters to their respective reactive zones
    !> The kinetic reaction rates of the initial state are zero (set at allocation): the interface
    !> routines must not rely on previous-step rates (see interfaz_comps_arch).
end subroutine