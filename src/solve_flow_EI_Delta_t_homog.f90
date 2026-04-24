subroutine solve_flow_EI_Delta_t_homog(this,theta,Time_out,output)
    !> Solves 1D transient flow equation with homogeneous time step using Euler implicit method
    
    !> this: transient PDE object
    !> theta: time weighting factor
    !> Time_out: output time values
    !> output: head vs time output
    
    !> Results at all intermediate steps are written in binary mode in file head_binary_EE.txt
    use flow_transient_m, only: flow_transient_c
    use time_discr_m, only: time_discr_homog_c
    use metodos_sist_lin_m, only: Thomas
    !use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs
    implicit none
    
    !> Variables
    class(flow_transient_c) :: this
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(out) :: output(:,:)

    integer(kind=4) :: n,i,icol,k,out_freq,head_r_flag,source_term_flag,Num_output
    real(kind=8) :: Time
    real(kind=8), parameter :: tol_Thomas=1d-9,epsilon=1d-9
    real(kind=8), allocatable :: head_old(:),head_new(:),b(:)
    !type(tridiag_matrix_c) :: E_mat,X_mat,A

    !procedure(Dirichlet_BCs_PDE), pointer :: p_BCs=>null()
    
    n=this%spatial_discr%Num_targets
    !select type (this)
    !class is (diffusion_1D_transient_c)
        select type (time_discr=>this%time_discr)
        type is (time_discr_homog_c)
            head_old=this%head_init
            allocate(head_new(n),b(n))
        !> BCs pointer
            !if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then
            !    !call Dirichlet_BCs_PDE(this)
            !    !p_BCs=>Dirichlet_BCs_PDE
            !else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then
            !    !call Neumann_homog_BCs(this)
            !    !p_BCs=>Neumann_homog_BCs
            !else
            !    error stop "Boundary conditions not implemented yet"
            !end if
            Num_output=size(Time_out)
        !> Implicit Euler
            open(unit=0,file="head_binary_EI.txt",form="unformatted",access="sequential",status="unknown") 
            icol=1
            Time=0
            if (abs(Time-Time_out(icol))<epsilon) then
                output(:,icol)=head_old
                icol=icol+1
            end if
            do k=1,time_discr%Num_time
                Time=k*time_discr%Delta_t
                write(0) Time, head_old
            !> Linear system
                call this%compute_b_vec_1D(theta,head_old,b,k)
                call Thomas(this%A_mat,b,tol_Thomas,head_new)
                if (abs(Time-Time_out(icol))<epsilon) then
                    output(:,icol)=head_new
                    icol=icol+1
                    if (icol>Num_output) then
                        write(*,*) "Reached Num_output"
                        exit
                    end if
                end if
                head_old=head_new
            end do
            this%head=head_new
            deallocate(head_old,head_new,b)
            close(0)
        end select
    !end select
end subroutine 