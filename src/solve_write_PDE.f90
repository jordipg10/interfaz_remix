subroutine solve_write_PDE(this,Time_out)
!> Calls subroutines that perform computations and writes results
    use PDE_m, only: PDE_c
    implicit none
    !> Variables
    class(PDE_c) :: this !> PDE object
    real(kind=8), intent(in) :: Time_out(:)
    
    integer(kind=4) :: Num_output
    real(kind=8), allocatable :: conc_out(:,:)
    
    Num_output=size(Time_out)
    allocate(conc_out(this%spatial_discr%Num_targets,Num_output))
    call this%solve_PDE()
    call this%write_PDE()
end subroutine
    
    


