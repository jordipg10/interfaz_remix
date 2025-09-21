subroutine solve_and_write_PDE_1D(this,root,Time_out)
!> Calls subroutines that perform computations and writes results
    use transport_transient_m
    use transport_m
    implicit none
    !> Variables
    class(PDE_1D_c) :: this !> PDE object
    character(len=*), intent(in) :: root
    real(kind=8), intent(in) :: Time_out(:)
    
    character(len=100) :: filename
    integer(kind=4) :: Num_time,Num_output, i,k,icol,out_freq
    real(kind=8) :: theta
    real(kind=8), allocatable :: conc_out(:,:)
    
    Num_output=size(Time_out)
    allocate(conc_out(this%spatial_discr%Num_targets,Num_output))
    call this%solve_PDE_1D(Time_out,conc_out)
    call this%write_PDE_1D(root,Time_out,conc_out)
end subroutine solve_and_write_PDE_1D
    
    


