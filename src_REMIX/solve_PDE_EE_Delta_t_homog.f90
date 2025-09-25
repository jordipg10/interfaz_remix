! subroutine solve_PDE_EE_Delta_t_homog(this,Time_out,output)
!     !> Solves 1D transient PDE with homogeneous time step using Lagr explicit method 
    
!     !> this: transient PDE object
!     !> Time_out: output time values
!     !> output: concentration vs time output
    
!     !> Results at all intermediate steps are written in binary mode in file conc_binary_EE.txt
    
!     use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
!     use transport_transient_m, only: transport_1D_transient_c, time_discr_homog_c, PDE_1D_transient_c, tridiag_matrix_c, &
!     mass_balance_error_ADE_trans_Dirichlet_recharge, mass_balance_error_ADE_trans_Dirichlet_discharge, &
!     mass_balance_error_ADE_trans_Dirichlet_evap, mass_balance_error_ADE_trans_PMF_recharge, &
!     mass_balance_error_ADE_trans_PMF_discharge, mass_balance_error_ADE_trans_PMF_evap
!     implicit none
    
!     !> Variables
!     class(PDE_1D_transient_c) :: this
!     real(kind=8), intent(in) :: Time_out(:)
!     real(kind=8), intent(out) :: output(:,:)

!     integer(kind=4) :: n,i,icol,k,out_freq,conc_r_flag,source_term_flag,Num_output
!     real(kind=8) :: Time,q_inf,q_out,mass_bal,mass_bal_err,cum_mass_bal_err,cum_sq_mass_bal_err,Delta_x
!     real(kind=8), parameter :: epsilon=1d-12
!     real(kind=8), allocatable :: conc_old(:),conc_new(:)
!     type(tridiag_matrix_c) :: E_mat,X_mat

!     !procedure(Dirichlet_BCs_PDE), pointer :: p_BCs=>null()
!     procedure(mass_balance_error_ADE_trans_Dirichlet_recharge), pointer :: p_mass_bal_err=>null()    
    
!     n=this%spatial_discr%Num_targets

!     select type (this)
!     type is (transport_1D_transient_c)
!         select type (time_discr=>this%time_discr)
!         type is (time_discr_homog_c)
!             conc_old=this%conc_init
!             allocate(conc_new(n))
!         !> BCs
!             if (this%BCs%BCs_label(1) == 1 .and. this%BCs%BCs_label(2) == 1) then
!                 !call Dirichlet_BCs_PDE(this)
!                 !p_BCs=>Dirichlet_BCs_PDE
!                 if (minval(this%tpt_props_heterog%source_term)>=0d0) then !> recharge
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_Dirichlet_recharge
!                 else if (maxval(this%tpt_props_heterog%source_term)<=0d0 .and.  this%BCs%evap.eqv..false.) then !> discharge
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_Dirichlet_discharge
!                 else !> evaporation
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_Dirichlet_evap
!                 end if
!             else if (this%BCs%BCs_label(1) == 2 .and. this%BCs%BCs_label(2) == 2) then
!                 !call Neumann_homog_BCs(this)
!             else if (this%BCs%BCs_label(1) == 3 .and. this%BCs%BCs_label(2) == 2) then
!                 !call Robin_Neumann_homog_BCs(this)
!                 !p_BCs=>Robin_Neumann_homog_BCs
!                 if (minval(this%tpt_props_heterog%source_term)>=0d0) then
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_PMF_recharge
!                 else if (maxval(this%tpt_props_heterog%source_term)<=0d0 .and.  this%BCs%evap.eqv..false.) then
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_PMF_discharge
!                 else
!                     p_mass_bal_err=>mass_balance_error_ADE_trans_PMF_evap
!                 end if
!             else
!                 error stop "Boundary conditions not implemented yet"
!             end if
!         !> Mixing ratios
!             !call this%compute_mixing_ratios(0d0)
!             !call this%compute_X_mat(0d0,X_mat)
!             !print *, X_mat%compute_norm_inf()
!         !> Explicit Euler
!             open(unit=0,file="conc_binary_EE.txt",form="unformatted",access="sequential",status="unknown")
!             Num_output=size(Time_out)
!             icol=1
!             Time=0
!             if (abs(Time-Time_out(icol))<epsilon) then
!                 output(:,icol)=conc_old
!                 icol=icol+1
!             end if
!             cum_mass_bal_err=0d0
!             Delta_x=this%spatial_discr%get_mesh_size()
!         !> Time integration loop
!             do k=1,time_discr%Num_time
!                 Time=k*time_discr%Delta_t
!                 write(0) Time, conc_old
!             !> Spatial discretization
!                 call this%compute_b_vec_lin_syst(0d0,conc_old,conc_new)
!             !> Mass balance error
!                 !mass_bal_err=p_mass_bal_err(this,conc_old,conc_new,time_discr%Delta_t,Delta_x)
!                 cum_mass_bal_err=cum_mass_bal_err + mass_bal_err
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
!             !print *, this%conc
!             deallocate(conc_old,conc_new)
!             close(0)
!         end select
!     end select
! end subroutine