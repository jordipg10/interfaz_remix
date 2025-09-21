! subroutine solve_PDE_EE_Delta_t_heterog(this,Time_out,output)
!     !> Solves 1D transient PDE with heterogeneous time step using Lagr explicit method
    
!     !> this: transient PDE object
!     !> Time_out: output time values
!     !> output: concentration vs time output
    
!     !> Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
!     use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
!     use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c, time_discr_heterog_c, PDE_1D_transient_c, &
!     tridiag_matrix_c
!     implicit none
    
!     !> Variables
!     class(PDE_1D_transient_c) :: this
!     real(kind=8), intent(in) :: Time_out(:)
!     real(kind=8), intent(out) :: output(:,:)

!     integer(kind=4) :: i,icol,k,out_freq,Num_output
!     real(kind=8) :: Time
!     real(kind=8), parameter :: epsilon=1d-9
!     real(kind=8), allocatable :: conc_old(:),conc_new(:)
    
!     type(tridiag_matrix_c) :: E_mat,X_mat

!     procedure(Dirichlet_BCs_PDE), pointer :: p_BCs=>null()
    
        
!     select type (this)
!     class is (diffusion_1D_transient_c)
!         select type (time_discr=>this%time_discr)
!         type is (time_discr_heterog_c)
!             conc_old=this%conc_init
!             allocate(conc_new(this%spatial_discr%Num_targets))
!         !> BCs pointer
!             if (this%BCs%BCs_label(1) == 1 .and. this%BCs%BCs_label(2) == 1) then
!                 !call Dirichlet_BCs_PDE(this)
!             else if (this%BCs%BCs_label(1) == 2 .and. this%BCs%BCs_label(2) == 2) then
!                 !call Neumann_homog_BCs(this)
!             else if (this%BCs%BCs_label(1) == 3 .and. this%BCs%BCs_label(2) == 2) then
!                 !call Robin_Neumann_homog_BCs(this)
!                 error stop "Boundary conditions not implemented yet"
!             end if
!         !> Explicit Euler
!             open(unit=0,file="conc_binary_EE.txt",form="unformatted",access="sequential",status="unknown")  
!             icol=1
!             Time=0
!             Num_output=size(Time_out)
!             if (abs(Time-Time_out(icol))<epsilon) then
!                 output(:,icol)=conc_old
!                 icol=icol+1
!             end if
!             do k=1,time_discr%Num_time
!                 Time=Time+time_discr%Delta_t(k)
!                 write(0) Time, conc_old
!             !> Spatial discretisation
!                 call this%compute_b_vec_lin_syst(0d0,conc_old,conc_new,k)
!                 if (abs(Time-Time_out(icol))<epsilon) then
!                     output(:,icol)=conc_new
!                     icol=icol+1
!                     if (icol>Num_output) then
!                         write(*,*) "Reached Num_output"
!                         exit
!                     end if
!                 end if
!                 conc_old=conc_new
!             end do
!             this%conc=conc_new
!             deallocate(conc_old,conc_new)
!             close(0)
!         end select
!     end select
! end subroutine 