subroutine main_PDE(this,root)
    use PDE_m, only: PDE_1D_c
    implicit none

!> Variables
    class(PDE_1D_c) :: this
    character(len=*), intent(in) :: root
    
    integer(kind=4) :: N_t
    real(kind=8), allocatable :: Time_out(:),Delta_t(:),Delta_r(:),Delta_r_D(:)
    real(kind=8) :: Final_time, measure, Delta_x, length, L2_norm_vi, radius
    real(kind=8) :: a, Delta_r_0
    integer(kind=4) :: scheme,int_method,Num_time,BCs,opcion_BCs,parameters_flag,i,j,Num_targets,eqn_flag,dim,niter_max,niter,&
        Dirichlet_BC_location,targets_flag
!****************************************************************************************************************************************************
!> Pre-process
    call this%initialise_PDE(root)
!****************************************************************************************************************************************************
!> Process
    !select type (this)
    !class is (diffusion_1D_transient_c)
    !    ! open(unit=5,file='C:\Users\user2319\source\repos\jordipg10\1D-Transport-Code\input\time_out.dat',status='old',action='read')
    !    ! read(5,*) N_t !> number of output times
    !    ! allocate(Time_out(N_t))
    !    ! read(5,*) Time_out !> output times
    !    ! close(5)
    !end select
    call this%solve_write_PDE_1D(Time_out)
!****************************************************************************************************************************************************
end subroutine