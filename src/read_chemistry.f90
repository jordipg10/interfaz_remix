!> This subroutine reads chemistry data from input files
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
    integer(kind=4) :: unit_chem_opts !> chemical options file unit
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    
    !> File units (arbitrary)
    unit_chem_syst=1
    unit_chem_opts=51
    unit_loc_chem=2

    !> Reads chemistry input files
    call this%read_chemistry_CHEPROO(root,path_pb,path_DB,unit_chem_syst,unit_loc_chem)
end subroutine