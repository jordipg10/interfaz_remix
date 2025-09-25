!> Computes variable activity species concentrations from component concentrations using Newton method
!! We assume all primary species are aqueous
!! We assume the initial guess for the primary concentrations is already set in the aqueous chemistry object
subroutine compute_c_nc_from_u_aq_Newton_ideal(this,conc_comp,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: conc_comp(:) !> component concentrations (molarities)
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations in Newton algorithm
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    real(kind=8), allocatable :: c1(:),c2nc(:),c2k(:),log_c2k(:),log_c2(:), dc2nc_dc1(:,:)
    real(kind=8), allocatable :: residual(:) !> c1+U2*c2nc-u
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: mat_lin_syst(:,:)
    real(kind=8), allocatable :: out_prod(:,:),d_log_gamma_d_I(:)
    integer(kind=4) :: i,n_e,n_p,niter_Picard,n_nc,n_nc_aq
    logical :: CV_flag_Picard,zero_flag
    
!> Pre-Process
    CV_flag=.false.
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_nc_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
    n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        
    allocate(c2nc(n_e),residual(n_p),dc2nc_dc1(n_e,n_p),Delta_c1(n_p),d_log_gamma_d_I(n_nc))
    niter=0
!> Process
        do
            niter=niter+1 !> we update number of iterations
            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
                print *, "Residual: ", inf_norm_vec_real(residual)
                print *, "Too many Newton iterations in speciation"
                exit
            end if
            call this%compute_c2nc_from_c1_aq_ideal(c2nc)
            conc_nc(1:n_p)=this%concentrations(this%ind_var_act_species(1:n_p)) !> chapuza
            conc_nc(n_p+1:n_nc)=c2nc !> chapuza
            call this%compute_residual(conc_comp,conc_nc,residual)
            if (inf_norm_vec_real(residual)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> CV reached
                CV_flag=.true.
                exit
            end if
        !> We compute Jacobian secondary-primary concentrations
            call this%compute_dc2nc_dc1_aq_ideal(c2nc,dc2nc_dc1)
        !> We solve linear system
            mat_lin_syst=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,1:n_p)+&
                matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,n_p+1:n_nc),dc2nc_dc1) !> U_1 + U_2_nc*dc2nc_dc1
            call LU_lin_syst(mat_lin_syst,-residual,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1)
            !> c1^(i+1)=c1^i+Delta_c1^i
            if (inf_norm_vec_real(Delta_c1/conc_nc(1:n_p))<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then !> chapuza
                print *, "Relative error: ", inf_norm_vec_real(Delta_c1/conc_nc(1:n_p))
                print *, "Newton speciation not accurate enough"
                exit
            else
                call this%update_conc_aq_prim_species(Delta_c1,zero_flag)
            end if
            !call this%compute_salinity()
            !call this%compute_molarities()
        end do
        !call this%compute_pH()
end subroutine