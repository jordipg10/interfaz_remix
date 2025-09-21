subroutine solve_tpt_1D_stat(this)
!> Solves 1D steady-state ADE using Thomas algorithm
    
    !> this: 1D tramsport object
        
    use transport_m, only: transport_1D_c, PDE_1D_c, diffusion_1D_c, mass_balance_error_ADE_stat_Dirichlet_discharge
    use metodos_sist_lin_m, only: tridiag_matrix_c, Thomas
    implicit none
    
    !> Variables
    class(transport_1D_c) :: this

    integer(kind=4) :: n,i,icol,k,out_freq,conc_r_flag,source_term_flag,Num_output
    real(kind=8) :: Time,MBE
    real(kind=8), parameter :: tol_Thomas=1d-9
    real(kind=8), allocatable :: conc_old(:),conc_new(:)
    type(tridiag_matrix_c) :: E_mat,X_mat

    procedure(mass_balance_error_ADE_stat_Dirichlet_discharge), pointer :: p_MBE=>null()
    
    n=this%spatial_discr%Num_targets
    
    !select type (this)
    !type is (transport_1D_c)
            allocate(this%conc(n))
            call Thomas(this%trans_mat,-this%source_term_PDE,tol_Thomas,this%conc) !> we solve linear system Tc=-g
            if (maxval(this%tpt_props_heterog%source_term)<0d0 .and. this%BCs%evap.eqv..false.) then !> discharge
                p_MBE=>mass_balance_error_ADE_stat_Dirichlet_discharge
            else
                error stop "Mass balance error not implemented yet"
            end if
            MBE=abs(p_MBE(this))
    !end select
end subroutine