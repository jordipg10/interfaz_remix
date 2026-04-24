subroutine write_flow_transient(this,root)
    use flow_transient_m, only: flow_transient_c
    implicit none
    class(flow_transient_c) :: this
    character(len=*), intent(in) :: root

    !> Write the flow transient object to a file
    !> This subroutine writes the flow transient object to a file for later use

    ! Open the file for writing
    open(unit=10, file=trim(root)//'_flow_transient.out', status='replace')

    ! Write the properties of the flow transient object
    ! write(10, *) 'Flow Transient Properties:'
    ! write(10, *) 'Time Step:', this%time_step
    ! write(10, *) 'Total Time:', this%total_time
    ! write(10, *) 'Number of Steps:', this%num_steps

    ! Close the file
    close(10)

end subroutine write_flow_transient