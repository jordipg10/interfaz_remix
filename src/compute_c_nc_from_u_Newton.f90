!> \file compute_c_nc_from_u_Newton.f90
!> \brief Computes variable activity species concentrations from component concentrations using Newton-Raphson method with non-ideal activity corrections
!> \details This subroutine solves the nonlinear speciation problem for non-ideal solutions (aqueous and solid)
!> using Newton-Raphson iteration with activity coefficient corrections. Given component concentrations and initial 
!> guesses for primary and secondary species, it computes both aqueous and solid secondary species concentrations 
!> by iteratively solving:
!> \f[
!> \mathbf{r}(\mathbf{c}_1) = \mathbf{c}_1 + \mathbf{U}_2 \cdot \mathbf{c}_{2nc} - \mathbf{u} = \mathbf{0}
!> \f]
!> where:
!> - \f$ \mathbf{c}_1 \f$ = primary species concentrations (aqueous + solid)
!> - \f$ \mathbf{c}_{2nc} \f$ = secondary non-component species concentrations
!> - \f$ \mathbf{U}_2 \f$ = component matrix for secondary species
!> - \f$ \mathbf{u} \f$ = component concentrations (input)
!>
!> Assumptions:
!> - Non-ideal solution behavior (activity coefficients γ ≠ 1)
!> - Initial guesses for both primary and secondary concentrations are provided
!> - Uses Picard iteration for c2v from c1 at each Newton step
!>
!> Newton iteration with activity corrections:
!> \f[
!> \mathbf{c}_1^{(i+1)} = \mathbf{c}_1^{(i)} + \Delta\mathbf{c}_1^{(i)}
!> \f]
!> where \f$ \Delta\mathbf{c}_1 \f$ solves:
!> \f[
!> \left(\mathbf{U}_1 + \mathbf{U}_2 \frac{\partial \mathbf{c}_{2nc}}{\partial \mathbf{c}_1}\right) \Delta\mathbf{c}_1 = -\mathbf{r}
!> \f]
!> with Jacobian including activity coefficient effects:
!> \f[
!> \frac{\partial \mathbf{c}_{2nc}}{\partial \mathbf{c}_1} = f(\mathbf{c}_1, \gamma, \frac{\partial \log \gamma}{\partial I})
!> \f]
!>
!> \param[in,out] this Aqueous chemistry object containing speciation algebra, solid chemistry, and concentrations
!> \param[in] c1_ig Initial guess for primary species concentrations (aqueous + solid) [mol/L_water]
!> \param[in] c2v_ig Initial guess for secondary non-component species concentrations [mol/L_water]
!> \param[in] conc_comp Component concentrations (total element concentrations) [mol/L_water]
!> \param[out] conc_nc Variable activity species concentrations: primary + secondary (already allocated) [mol/L_water]
!> \param[out] niter Number of Newton iterations performed [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE otherwise [-]

subroutine compute_c_nc_from_u_Newton(this,c1_ig,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real, outer_prod_vec !> Import aqueous chemistry class and utility functions for infinity norm, linear system solver, and outer product
    use metodos_sist_lin_m, only: LU_lin_syst !> Import linear system solver using LU decomposition
    implicit none !> Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone, speciation algebra, and species concentrations (aqueous + solid) [-]
    real(kind=8), intent(in) :: c1_ig(:) !> Initial guess for primary species concentrations (used to initialize Newton iteration) [mol/L_water]
    real(kind=8), intent(in) :: c2v_ig(:) !> Initial guess for secondary variable activity species concentrations (used for Picard iteration) [mol/L_water]
    real(kind=8), intent(in) :: conc_comp(:) !> Component concentrations (total element concentrations) [mol/L_water]
    real(kind=8), intent(out) :: conc_nc(:) !> Variable activity species concentrations: primary + secondary (already allocated) [mol/L_water]
    integer(kind=4), intent(out) :: niter !> Number of iterations performed in Newton-Raphson algorithm [-]
    logical, intent(out) :: CV_flag !> Convergence flag: TRUE if converged within tolerance, FALSE if failed to converge [-]
!> Variables
    real(kind=8), allocatable :: c1(:),c2v_old(:),c2v_new(:),c2k(:),log_c2k(:),log_c2(:), dc2v_dc1(:,:) !> c1: primary species concentrations (aqueous + solid) [mol/L_water], c2v_old: previous secondary concentrations for Picard iteration [mol/L_water], c2v_new: new secondary concentrations [mol/L_water], c2k: kinetic species concentrations (unused) [mol/L_water], log_c2k: log kinetic concentrations (unused) [-], log_c2: log secondary concentrations (unused) [-], dc2v_dc1: Jacobian matrix ∂c2v/∂c1 including activity coefficient effects [mol/mol]
    real(kind=8), allocatable :: residual(:) !> Residual vector r = c1 + U2·c2v - u (component mass balance error) [mol/L_water]
    real(kind=8), allocatable :: Delta_c1(:) !> Newton correction: Δc1 = c1^(i+1) - c1^(i) (change in primary concentrations) [mol/L_water]
    real(kind=8), allocatable :: mat_lin_syst(:,:) !> Linear system matrix: U1 + U2·(∂c2v/∂c1) for Newton iteration including activity corrections [-]
    real(kind=8), allocatable :: out_prod(:,:),out_prod_aq(:,:),d_log_gamma_d_I(:),log_Jacobian_act_coeffs(:,:),& !> out_prod: outer product matrix of ∂(log γ)/∂I and z² for all species [-], out_prod_aq: aqueous-only outer product [-], d_log_gamma_d_I: derivative of log(activity coefficient) w.r.t. ionic strength [-/M], log_Jacobian_act_coeffs: log-Jacobian of activity coefficients w.r.t. concentrations [-]
        log_Jacobian_act_coeffs_aq(:,:) !> log_Jacobian_act_coeffs_aq: log-Jacobian of aqueous activity coefficients w.r.t. aqueous concentrations [-]
    integer(kind=4) :: i,n_e,n_p,n_p_aq,niter_Picard,n_v,n_v_aq !> i: loop index [-], n_e: number of equilibrium reactions [-], n_p: number of primary species (aqueous + solid) [-], n_p_aq: number of aqueous primary species [-], niter_Picard: number of Picard iterations for c2v computation [-], n_v: total number of variable activity species [-], n_v_aq: number of aqueous variable activity species [-]
    logical :: CV_flag_Picard !> CV_flag_Picard: convergence flag for Picard iteration (TRUE if Picard converged) [-]
    
!> Pre-Process: Initialize variables and allocate arrays
    CV_flag=.false. !> Initialize convergence flag to FALSE (assume not converged) [-]
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !> Extract number of equilibrium reactions from speciation algebra [-]
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Extract number of primary species (includes aqueous and solid) from speciation algebra [-]
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> Extract number of aqueous primary species only from speciation algebra [-]
    n_v_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !> Extract number of aqueous variable activity species from speciation algebra [-]
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !> Extract total number of variable activity species (aqueous + solid) from speciation algebra [-]
        
    allocate(c2v_new(n_e),residual(n_p),dc2v_dc1(n_e,n_p),Delta_c1(n_p),d_log_gamma_d_I(n_v),out_prod_aq(n_v_aq,n_v_aq),& !> Allocate arrays: c2v_new(n_e) for new secondary concentrations, residual(n_p) for component balance error, dc2v_dc1(n_e,n_p) for Jacobian, Delta_c1(n_p) for Newton correction, d_log_gamma_d_I(n_v) for activity coefficient derivatives, out_prod_aq(n_v_aq,n_v_aq) for aqueous outer products
    log_Jacobian_act_coeffs_aq(n_v_aq,n_v_aq),log_Jacobian_act_coeffs(n_v,n_v)) !> log_Jacobian_act_coeffs_aq(n_v_aq,n_v_aq) for aqueous log-Jacobian, log_Jacobian_act_coeffs(n_v,n_v) for full log-Jacobian
    niter=0 !> Initialize Newton iteration counter to zero [-]
    log_Jacobian_act_coeffs=0d0 !> Initialize log-Jacobian of activity coefficients to zero (temporary workaround) [-]
    d_log_gamma_d_I=0d0 !> Initialize derivative of log(activity coefficient) w.r.t. ionic strength to zero (temporary workaround) [-/M]
    conc_nc(1:n_p)=c1_ig !> Initialize primary species concentrations in output array with initial guess [mol/L_water]
    c2v_old=c2v_ig !> Store initial guess for secondary concentrations to use in first Picard iteration (temporary workaround) [mol/L_water]
!> Process: Newton-Raphson iteration loop with activity coefficient corrections
        do !> Begin infinite loop (will exit when converged or max iterations reached)
            niter=niter+1 !> Increment Newton iteration counter [-]
            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then !> Check if maximum number of iterations exceeded
                print *, "Residual: ", inf_norm_vec_real(residual) !> Print final residual infinity norm for diagnostics [mol/L_water]
                print *, "Too many Newton iterations in speciation" !> Print error message indicating failure to converge
                exit !> Exit iteration loop due to exceeding maximum iterations
            end if
            call this%compute_c2v_from_c1_Picard(conc_nc(1:n_p),c2v_old,conc_nc(n_p+1:n_v),niter_Picard,CV_flag_Picard) !> Compute secondary non-component concentrations from current primary concentrations using Picard iteration with activity coefficient corrections (uses c2v_old as initial guess) [mol/L_water]
            call this%compute_res_spec(conc_comp,conc_nc,residual) !> Compute residual vector: r = c1 + U2·c2v - u (component mass balance error) [mol/L_water]
            if (inf_norm_vec_real(residual)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> Check if infinity norm of residual is below absolute tolerance (convergence criterion)
                CV_flag=.true. !> Set convergence flag to TRUE (solution found) [-]
                exit !> Exit iteration loop due to successful convergence
            end if
        !> First we compute d_log_gamma_d_I
            !d_log_gamma_d_I=this%compute_d_log_gamma_d_I_aq_chem()
        !> Outer product d_log_gamma_nc_d_I and z_nc^2
            !out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:n_v))
        !> We compute Jacobian secondary variable activity-primary concentrations
            call this%compute_dc2v_dc1(conc_nc(1:n_p),conc_nc(n_p+1:n_v),dc2v_dc1)
        !> We compute log-Jacobian variable activity coefficients-variable activity concentrations
            !out_prod_aq(1:n_p_aq,1:n_p_aq)=out_prod(1:n_p_aq,1:n_p_aq)
            !out_prod_aq(n_p_aq+1:n_v_aq,n_p_aq+1:n_v_aq)=out_prod(n_p+1:n_v_aq+1,n_p+1:n_v_aq+1)
            !call this%compute_log_Jacobian_act_coeffs_aq_chem()            
        !> We check Jacobain secondary variable activity-primary concentrations
             !call this%check_dc2v_dc1(conc_nc(1:n_p),conc_nc(n_p+1:n_v),dc2v_dc1,log_Jacobian_act_coeffs)
        !> We solve linear system
            mat_lin_syst=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,1:n_p)+matmul(&
                this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,n_p+1:n_v),dc2v_dc1) !> U_1 + U_2_nc*dc2v_dc1
            call LU_lin_syst(mat_lin_syst,-residual,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1)
            !> c1^(i+1)=c1^i+Delta_c1^i
            if (inf_norm_vec_real(Delta_c1/conc_nc(1:n_p))<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then !> chapuza
                print *, "Relative error: ", inf_norm_vec_real(Delta_c1/conc_nc(1:n_p))
                print *, "Newton speciation not accurate enough"
                exit
            else
                call this%update_conc_prim_species(conc_nc(1:n_p),Delta_c1)
                c2v_old=conc_nc(n_p+1:n_v)
            end if
            !call this%compute_salinity()
            !call this%compute_molarities()
        end do
        call this%compute_pH()
end subroutine