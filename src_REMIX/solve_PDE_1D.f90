! !> Calls subroutines that perform 1D PDE computations
! subroutine solve_PDE_1D(this,Time_out,output)
!     use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs
!     use PDE_transient_m, only: PDE_1D_transient_c, PDE_1D_c
!     use diffusion_transient_m, only: diffusion_1D_transient_c
!     use diffusion_m, only: diffusion_1D_c
!     use time_discr_m, only: time_discr_homog_c, time_discr_c, time_discr_heterog_c
!     use matrices_m, only: tridiag_matrix_c  
!     implicit none
! !> Variables
!     class(PDE_1D_c) :: this !> 1D PDE object
!     real(kind=8), intent(in) :: Time_out(:)
!     real(kind=8), intent(out) :: output(:,:)
    
!     integer(kind=4) :: k,n,j,i
!     real(kind=8) :: sum,B_norm_inf
!     real(kind=8), allocatable :: Delta_r(:)
!     real(kind=8), parameter :: tol=1d-12
!     type(tridiag_matrix_c) :: X_mat,A_mat
    
!     n=this%spatial_discr%Num_targets
    
! !> We compute arrays
!     call this%allocate_arrays_PDE_1D()
!     select type (this)
!     class is (diffusion_1D_c)
!         if (this%sol_method.eq.1) then
!             call this%compute_trans_mat_PDE()
!             call this%compute_source_term_PDE()
!             call this%solve_PDE_1D_stat()
!         end if
!     class is (PDE_1D_transient_c)
!         if (this%sol_method.eq.1) then
!             select type (time_discr=>this%time_discr)
!             type is (time_discr_homog_c)
!                 call this%compute_mixing_ratios_Delta_t_homog()
!                 if (time_discr%int_method.eq.1) then
!                     !call this%compute_mixing_ratios_Delta_t_homog(0d0)
!                     call this%solve_PDE_EE_Delta_t_homog(Time_out,output)
!                 else if (time_discr%int_method.eq.2) then
!                     !call this%compute_mixing_ratios_Delta_t_homog(1d0)
!                     call this%solve_PDE_EI_Delta_t_homog(1d0,Time_out,output)
!                 else if (time_discr%int_method.eq.3) then
!                     !call this%compute_mixing_ratios_Delta_t_homog(5d-1)
!                     call this%solve_PDE_EI_Delta_t_homog(5d-1,Time_out,output)
!                 else if (time_discr%int_method.eq.4) then
!                     call this%solve_PDE_RKF45(time_discr%Delta_t,tol)
!                 end if
!             type is (time_discr_heterog_c)
!                 if (time_discr%int_method.eq.1) then
!                     call this%solve_PDE_EE_Delta_t_heterog(Time_out,output)
!                 end if
!             end select
!         else if (this%sol_method.eq.2) then
!             select type (this)
!             type is (diffusion_1D_transient_c)
!                 call this%compute_A_mat_ODE(A_mat)
!                 call A_mat%compute_eigenvalues()
!                 call A_mat%compute_eigenvectors()
!                 if (minval(A_mat%eigenvalues)<=0d0) then
!                     error stop "Eigenvalues are not positive"
!                 end if
!                 allocate(A_mat%eigenvectors(n,n))
!                 call A_mat%check_eigenvectors_tridiag_sym_matrix(tol)
!             end select
!         end if
!     end select
! end subroutine solve_PDE_1D