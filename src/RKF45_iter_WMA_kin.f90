!!> f(c)=L(c)*inv(F)+transp(S_k)*rk(c)*inv(F)=(c*transp(T)+c_ext+transp(S_k)*rk(c))*inv(F)
!subroutine RKF45_iter_WMA_kin(this,F_mat,trans_mat,Delta_t_old,conc_old,conc_ext,tolerance,conc_new,Delta_t_new,niter,k)
!>    use chemistry_kin_m
!>    !use matrices_m
!>    implicit none
!>    class(chemistry_kin_c) :: this
!>    !integer(kind=4), intent(in) :: j !> target
!>    !class(tridiag_matrix_vec_c), intent(in) :: mixing_ratios
!>    class(diag_matrix_c), intent(in) :: F_mat
!>    class(tridiag_matrix_c), intent(in) :: trans_mat !> T
!>    real(kind=8), intent(in) :: Delta_t_old
!>    real(kind=8), intent(in) :: conc_old(:,:) !> c^k
!>    real(kind=8), intent(in) :: conc_ext(:,:) !> c_e
!>    real(kind=8), intent(in) :: tolerance
!>    real(kind=8), intent(out) :: conc_new(:,:) !> c^(k+1)
!>    real(kind=8), intent(out) :: Delta_t_new
!>    integer(kind=4), intent(out) :: niter !> number of iterations !> number of iterations
!>    integer(kind=4), intent(in), optional :: k
!>    
!>    integer(kind=4) :: n,i,icol,Num_output,counter
!>    real(kind=8) :: Time,Delta_t_max,Delta_t
!>    real(kind=8), allocatable :: k_RKF(:,:,:),conc_RK4_old(:,:),conc_RK4_new(:,:),conc_RK5(:,:)
!>    
!>    n=size(this%aqueous_chem)
!>    !allocate(k(n,6))
!>    
!>    !select type (this)
!>    !type is (chemistry_transient_c)
!>        !select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>        !type is (chem_system_eq_kin_c)
!>            !do j=1,n
!>                !allocate(conc_RK5(n))
!>                conc_RK4_old=conc_old
!>                Delta_t=Delta_t_old
!>                Delta_t_max=0d0
!>                !Time=Delta_t_init
!>                !this%time_discr%Num_time=0
!>                !icol=1
!>                niter=0
!>                do
!>                    niter=niter+1 !> we update number of iterations
!>                    k_RKF=this%compute_k_RKF45_RT_kin(F_mat,trans_mat,Delta_t,conc_RK4_old,conc_ext)
!>                    conc_RK4_new=conc_RK4_old+Delta_t_old*(25*k_RKF(:,:,1)/216+1408*k_RKF(:,:,3)/2565+2197*k_RKF(:,:,4)/4104-k_RKF(:,:,5)/5)
!>                    conc_RK5=conc_RK4_old+Delta_t*(16*k_RKF(:,:,1)/135+6656*k_RKF(:,:,3)/12825+28561*k_RKF(:,:,4)/56430-9*k_RKF(:,:,5)/50+2*k_RKF(:,:,6)/55)
!>                    Delta_t_new=9d-1*Delta_t*(tolerance/norm_mat_inf(conc_RK4_new-conc_RK5))**0.2 !> update time step
!>                    !print *, Delta_t_new
!>                    if (Delta_t_new>Delta_t_max) then
!>                        Delta_t_max=Delta_t_new
!>                    end if
!>                    if (norm_mat_inf(conc_RK4_new-conc_RK5)<tolerance) then
!>                        !Time=Time+Delta_t_new
!>                        conc_RK4_old=conc_RK4_new
!>                        conc_new=conc_RK5
!>                        !conc_comp_new=matmul(chem_syst%mob_comp_mat,conc_RK5)
!>                        exit
!>                    end if
!>                    Delta_t=Delta_t_new
!>                end do
!>                !this%conc=conc_RK5
!>        
!>        
!>        !end select
!>    !print *, Delta_t_max
!end subroutine