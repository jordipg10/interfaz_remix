!> Lectura quimica
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

    !> Read chemical options
    call this%read_chem_opts(path_pb//root,unit_chem_opts) 
    !> Read chemistry depending on the selected option
    if (this%read_opt==1) then !> CHEPROO-based
        call this%read_chemistry_CHEPROO(root,path_pb,path_DB,unit_chem_syst,unit_loc_chem)
    else if (this%read_opt==2) then !> PHREEQC
        error stop "PHREEQC data input not fully implemented yet"
    else if (this%read_opt==3) then !> PFLOTRAN
        error stop "PFLOTRAN data input not fully implemented yet"
    else
        error stop "This data input option is not implemented yet"
    end if
end subroutine