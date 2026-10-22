!> \file read_chemistry.f90
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
subroutine read_chemistry(this,root,path_pb,path_DB)
    use chemistry_m, only: chemistry_c
    implicit none
    !> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root of input and output files of problem
    character(len=*), intent(in) :: path_pb !> path for input and output files of problem
    character(len=*), intent(in) :: path_DB !> path for chemical databases

    !> Variables
    integer(kind=4) :: unit_chem_syst !> chemical system file unit
    integer(kind=4) :: unit_comp_opts !> computation options file unit
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    integer(kind=4) :: read_opt !> option for reading chemistry data
    integer(kind=4) :: i !> reactive zone loop index
    integer(kind=4), allocatable :: twiz_loc(:) !> domain target-water indices of one reactive zone
    
    !> File units (arbitrary)
    unit_chem_syst=1
    unit_comp_opts=51
    unit_loc_chem=2

    !> Read computation options
    call this%read_comp_opts(path_pb,root,unit_comp_opts)
    !> Read convergence parameters
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
    !> Assign each reactive zone its canonical id (position in reactive_zones) here, so that
    !> set_wat_indices_zones can match the waters to the zones by id. The grouping of the waters
    !> by reactive zone (wat_indices_zones) is NOT built here: set_target_waters (RT_m) may still
    !> rebuild the waters array (edge-centred Dirichlet case: the boundary waters are appended to
    !> the target waters and every water index shifts), so it is built at the end of
    !> set_target_waters, once the water indices are final.
    do i=1,this%num_reactive_zones
        this%reactive_zones(i)%id=i
    end do
    ! call this%wat_indices_zones%allocate_array(this%num_reactive_zones)
    ! do i=1,this%num_reactive_zones
    !     call this%set_wat_indices_zones(i,twiz_loc)
    !     call this%wat_indices_zones%cols(i)%allocate_vector(size(twiz_loc))
    !     this%wat_indices_zones%cols(i)%col_1=twiz_loc
    ! end do
    !> Read output options
     call this%chem_out_options%read_chem_out_options(path_pb,root,this%chem_syst,this%num_waters,this%num_rech_waters)
end subroutine