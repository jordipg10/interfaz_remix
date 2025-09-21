! !> Solves 1D transient PDE using RKF45 method
! subroutine solve_PDE_RKF45(this,Delta_t_init,tolerance)
!     use PDE_transient_m, only: PDE_1D_transient_c
!     use transport_transient_m, only: diffusion_1D_transient_c, time_discr_homog_c
!     use vectors_m
!     implicit none
!     class(PDE_1D_transient_c) :: this
!     real(kind=8), intent(in) :: Delta_t_init
!     real(kind=8), intent(in) :: tolerance
    
!     integer(kind=4) :: n,i,j,icol,Num_time,counter
!     real(kind=8) :: Delta_t_old,Delta_t_new,Time,Delta_t_max
!     real(kind=8), allocatable :: k(:,:),conc_RK4_old(:),conc_RK4_new(:),conc_RK5(:)
    
!     n=this%spatial_discr%Num_targets
    
!     select type (this)
!     class is (diffusion_1D_transient_c)
!         conc_RK4_old=this%conc_init
!         Delta_t_old=Delta_t_init
!         Delta_t_max=0d0 !> initialisation
!         Time=Delta_t_init
!         Num_time=0
!         icol=1
!         counter=0
!         do while(Time<this%time_discr%Final_time)
!             counter=counter+1
!             k=this%compute_k_RKF45(Delta_t_old,conc_RK4_old)
!             conc_RK4_new=conc_RK4_old+Delta_t_old*(25*k(:,1)/216+1408*k(:,3)/2565+2197*k(:,4)/4104-k(:,5)/5)
!             conc_RK5=conc_RK4_old+Delta_t_old*(16*k(:,1)/135+6656*k(:,3)/12825+28561*k(:,4)/56430-9*k(:,5)/50+2*k(:,6)/55)
!             Delta_t_new=9d-1*Delta_t_old*(tolerance/p_norm_vec(conc_RK4_new-conc_RK5,2))**0.2 !> update time step
!             if (Delta_t_new>Delta_t_max) then
!                 Delta_t_max=Delta_t_new
!             end if
!             if (p_norm_vec(conc_RK4_new-conc_RK5,2)<tolerance) then
!                 Num_time=Num_time+1
!                 Time=Time+Delta_t_new
!                 conc_RK4_old=conc_RK4_new
!             end if
!             Delta_t_old=Delta_t_new
!         end do
!         this%conc=conc_RK5
!         call this%time_discr%set_Num_time(Num_time)
!         call this%time_discr%set_Final_time(Time)
!         deallocate(k,conc_RK4_old,conc_RK4_new,conc_RK5)
!     end select
! end subroutine