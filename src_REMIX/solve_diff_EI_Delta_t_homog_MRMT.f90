subroutine solve_diff_EI_Delta_t_homog_MRMT(this,theta,Time_out,output)
    !> Solves 1D transient diffusion PDE with homogeneous time step using Euler implicit method and MRMT
    
    !> this: transient diffusion object
    !> theta: time weighting factor
    !> Time_out: output time values
    !> output: concentration vs time output
    
    !> Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
    use MRMT_m, only: MRMT_c, diffusion_1D_transient_c, time_discr_homog_c, tridiag_matrix_c, Thomas
    implicit none
    
    !> Variables
    class(MRMT_c) :: this
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(out) :: output(:,:)

    integer(kind=4) :: n,i,icol,k,out_freq,conc_r_flag,source_term_flag,Num_output,n_imm
    real(kind=8) :: Time,phi_mob,phi_imm,D,gk,sum_gk2,cst
    real(kind=8), parameter :: epsilon=1d-9,tol_Thomas=1d-9
    real(kind=8), allocatable :: conc_mob_old(:),conc_mob_new(:),conc_imm_init(:),conc_imm_old(:),conc_imm_new(:),b(:),alpha(:),&
    prob(:)
    type(tridiag_matrix_c) :: A_mat,A_mat_ODE
    
    n=this%PDE%spatial_discr%Num_targets

    !open(unit=56,file='dist_flujos.out',status='unknown')
    select type (PDE=>this%PDE)
    class is (diffusion_1D_transient_c)
        select type (time_discr=>PDE%time_discr)
        type is (time_discr_homog_c)
            conc_mob_old=PDE%conc_init
            conc_mob_old(n)=PDE%BCs%conc_out
            call this%set_n_imm(n) !> each target associated to one immobile zone
            this%mob_zone%mob_por=4d-1
            allocate(conc_mob_new(n),conc_imm_init(this%n_imm),conc_imm_old(this%n_imm),conc_imm_new(this%n_imm),b(n))
            conc_imm_init=0d0
            conc_imm_old=conc_imm_init
            call this%allocate_imm_zones()
            call PDE%compute_A_mat_ODE(A_mat_ODE)
            call A_mat_ODE%compute_eigenvalues()
            !cst=0d0
            !do i=1,n
            !>    cst=cst+1d0/(A_mat_ODE%eigenvalues(i)**pot)
            !end do
            !cst=1d0/cst
            !do i=1,this%n_imm
            !>    this%imm_zones(i)%exch_rate=A_mat_ODE%eigenvalues(i)
            !>    this%imm_zones(i)%prob=cst/(A_mat_ODE%eigenvalues(i)**pot)
            !>    print *, this%imm_zones(i)%prob
            !>    this%imm_zones(i)%imm_por=1d0
            !end do
            call this%check_imm_zones()
            call this%compute_A_mat_conc_mob(theta,time_discr%Delta_t,A_mat)
            Num_output=size(Time_out)
            open(unit=0,file="conc_binary_EI_MRMT.txt",form="unformatted",access="sequential",status="unknown")
            icol=1
            Time=0
            if (abs(Time-Time_out(icol))<epsilon) then
                output(:,icol)=[conc_mob_old,conc_imm_old]
                icol=icol+1
            end if
            !sum_gk2=0d0
            do k=1,time_discr%Num_time
                Time=k*time_discr%Delta_t
                write(0) Time, [conc_mob_old,conc_imm_old]
            !> Linear system
                call this%compute_b_conc_mob(theta,time_discr%Delta_t,conc_mob_old,conc_imm_old,b)
                call Thomas(A_mat,b,tol_Thomas,conc_mob_new)
                call this%compute_conc_imm_MRMT(theta,conc_imm_old,conc_mob_old,conc_mob_new,time_discr%Delta_t,conc_imm_new)
                if (abs(Time-Time_out(icol))<epsilon) then
                    output(:,icol)=[conc_mob_new,conc_imm_new]
                    icol=icol+1
                    if (icol>Num_output) then
                        write(*,*) "Reached Num_output"
                        exit
                    end if
                end if
                conc_mob_old=conc_mob_new
                conc_imm_old=conc_imm_new
                call PDE%prod_total_conc(A_mat_ODE,Time)
                !sum_gk2=sum_gk2+gk**2
            end do
            !write(56,*) sum_gk2
            this%conc_mob=conc_mob_new
            do i=1,this%n_imm
                this%conc_imm(i)=conc_imm_new(i)
            end do
            deallocate(conc_mob_old,conc_mob_new,conc_imm_old,conc_imm_new)
            close(0)
        end select
    end select
    !close(56)
end subroutine 