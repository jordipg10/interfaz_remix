!> Writes data and final results of 1D reactive transport problem
subroutine write_RT_1D(this,root,path_py)
    use RT_1D_m, only: RT_1D_c, RT_1D_transient_c
    implicit none
    class(RT_1D_c), intent(in) :: this                      !> 1D reactive transport object
    !integer(kind=4), intent(in) :: unit                     !> unit of output file
    character(len=*), intent(in) :: root                    !> root of output file
    character(len=*), intent(in), optional :: path_py       !> path for Python files
    
    integer(kind=4) :: unit
    
    unit=6 !> arbitrary unit for output file
    
    open(unit,file=root//'.out',status='unknown',form='formatted')
    select type (this)
    type is (RT_1D_transient_c)
    !> First we write transport data
        call this%transport%write_transport_data_WMA(unit)
    !> Then we write chemical data & results
        call this%chemistry%write_chemistry(unit,this%transport%Lagr_flag)
    !> Optionally, we write data for Python
        if (present(path_py)) then
            call this%write_python(path_py)
        end if
        !if (this%transport%Lagr_flag) then
        !    !> Write time discretisation data
        !    call this%write_Lagr_info(unit)
        !end if
    end select
    close(unit)
end subroutine write_RT_1D