!!> x=U_m*c_nc*F
!!> f=U_m*L(c_nc)+U_m*S_k,nc^T*r_k(c)
!subroutine RKF45_iter_WMA(this,mixing_ratios,F_mat,Delta_t_old,conc_old,conc_comp_new,Delta_t_new,niter,k)
!>    use chemistry_m
!>    use vectores_m
!>    implicit none
!>    class(chemistry_c) :: this
!>    !integer(kind=4), intent(in) :: j !> target
!>    class(tridiag_matrix_vec_c), intent(in) :: mixing_ratios
!>    class(diag_matrix_c), intent(in) :: F_mat
!>    real(kind=8), intent(in) :: Delta_t_old
!>    real(kind=8), intent(in) :: conc_old(:,:) !> (c^k|c_e)
!>    real(kind=8), intent(inout) :: conc_comp_new(:,:) !> u_m^(k+1)
!>    real(kind=8), intent(out) :: Delta_t_new
!>    integer(kind=4), intent(out) :: niter !> number of iterations !> number of iterations
!>    integer(kind=4), intent(in), optional :: k
!>    
!>    integer(kind=4) :: n,i,j,icol,Num_output,counter
!>    real(kind=8) :: Delta_t_old,Delta_t_new,Time,Delta_t_max
!>    real(kind=8), allocatable :: k(:,:),conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:)
!>    
!>    n=size(this%aqueous_chem)
!>    !allocate(k(n,6))
!>    
!>    select type (this)
!>    type is (chemistry_transient_c)
!>        select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>        type is (chem_system_eq_kin_c)
!>            !do j=1,n
!>                !allocate(conc_RK5(n))
!>                conc_RK4_old=conc_old(:,1:n)
!>                !Delta_t_old=Delta_t_init
!>                !Delta_t_max=0d0
!>                !Time=Delta_t_init
!>                !this%time_discr%Num_time=0
!>                !icol=1
!>                !niter=0
!>                !Num_output=size(Time_out)
!>                !if (abs(Time-Time_out(icol))<tolerance) then
!>                !>    output(:,icol)=conc_RK4_old
!>                !>    icol=icol+1
!>                !end if
!>                do !while(Time<this%time_discr%Final_time)
!>                    !if (abs(Time-Time_out(icol))<tolerance) then
!>                    !>    output(:,icol)=conc_RK5
!>                    !>    icol=icol+1
!>                    !>    if (icol>Num_output) then
!>                    !>        write(*,*) "Reached Num_output"
!>                    !>        exit
!>                    !>    end if
!>                    !end if
!>                    niter=niter+1 !> we update number of iterations
!>                    k=this%compute_k_RKF45(Delta_t_old,conc_RK4_old)
!>                    conc_RK4_new=conc_RK4_old+Delta_t_old*(25*k(:,1)/216+1408*k(:,3)/2565+2197*k(:,4)/4104-k(:,5)/5)
!>                    conc_RK5=conc_RK4_old+Delta_t_old*(16*k(:,1)/135+6656*k(:,3)/12825+28561*k(:,4)/56430-9*k(:,5)/50+2*k(:,6)/55)
!>                    Delta_t_new=9d-1*Delta_t_old*(tolerance/p_norm_vec(conc_RK4_new-conc_RK5,2))**0.2 !> update time step
!>                    !print *, Delta_t_new
!>                    if (Delta_t_new>Delta_t_max) then
!>                        Delta_t_max=Delta_t_new
!>                    end if
!>                    if (p_norm_vec(conc_RK4_new-conc_RK5,2)<tolerance) then
!>                        !this%time_discr%Num_time=this%time_discr%Num_time+1
!>                        !Time=Time+Delta_t_new
!>                        conc_RK4_old=conc_RK4_new
!>                        conc_comp_new=matmul(chem_syst%mob_comp_mat,conc_RK5)
!>                        exit
!>                    end if
!>                    Delta_t_old=Delta_t_new
!>            
!>            
!>                end do
!>                !this%conc=conc_RK5
!>        
!>        
!>        end select
!>    print *, Delta_t_max
!end subroutine