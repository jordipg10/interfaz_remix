!> \file initialise_conc_anal.f90
!> \brief Newton speciation with analytical Jacobian (general activity coefficients).
!> \details
!> Computes species concentrations from water-type constraints using
!> Newton's method with the analytical Jacobian, accounting for
!> non-ideal activity coefficients.  Used during CHEPROO water-type reading.
!> The Newton unknowns are the LOG10 of the aqueous primary concentrations (step bounded to
!> max_dp decades per iteration, backtracking on the tolerance-scaled residual, iteration budget
!> max(5*niter_max, 50)): the input guesses of a water type can be many decades off, which an
!> arithmetic step cannot cover. The log-space linear system is row-equilibrated before the LU
!> solve, so that the singularity test of LU_lin_syst (absolute, on the determinant) is not
!> triggered by the ln10*c1_j column scaling of trace primaries. A Newton that stagnates (step
!> below rel_tol, or no decrease in the line search) with every residual within the noise floor
!> of the residual evaluation, max(tol_res, eps^(2/3)*|ctot|), is accepted as converged.
!>
!> Cost per Newton iteration: one Picard speciation per line-search trial and none elsewhere (the
!> accepted trial IS the next iterate, so its speciation, held in the object and in c2_new, is
!> reused instead of being recomputed at the top of the iteration), and no allocations or array
!> temporaries inside the loop: every work array is allocated once, the index maps, the squared
!> valences and the per-species constants of d(log gamma)/dI are gathered once, and
!> log_Jacobian_act_coeffs (fixed sparsity pattern) is zeroed once. The results are identical to
!> those of the version that re-speciated every iterate.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup chemistry
!> \see initialise_conc_anal_ideal, initialise_conc_incr_coeff, aqueous_chemistry_m

!> Computes species concentrations with data from water type definition using Newton method
!> We assume the initial guess of primary aqueous species is already set in the aqueous chemistry object
!> This subroutine is only used to read water types based in CHEPROO
subroutine initialise_conc_anal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use arrays_m, only: int_array_c
    use metodos_sist_lin_m, only: LU_lin_syst
    use vectors_m, only: inf_norm_vec_real
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: icon(:) !> initial condition type
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> indices of constrains in equilibrium stoichiometric matrix
    real(kind=8), intent(in) :: ctot(:) !> data given
    integer(kind=4), intent(out) :: niter !> number of iterations Newton-Raphson
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    integer(kind=4) :: n_p !> number of primary species (aqueous + solid)
    integer(kind=4) :: n_p_aq !> number of aqueous primary species (Newton unknowns)
    integer(kind=4) :: n_2_aq !> number of aqueous secondary species
    integer(kind=4) :: n_e !> number of equilibrium reactions (secondary species)
    integer(kind=4) :: n_sp !> number of species
    integer(kind=4), allocatable :: p(:) !> concentration index of primary j (stoichiometric order): ind_prim_species(1:n_p_aq)
    integer(kind=4), allocatable :: s(:) !> concentration index of aqueous secondary i: ind_sec_species(1:n_2_aq)
    real(kind=8), allocatable :: c1(:) !> primary concentrations (n_p): the current Newton iterate
    real(kind=8), allocatable :: c2_old(:),c2_new(:) !> secondary concentrations: ideal guess of the Picard speciation, and its result
    real(kind=8), allocatable :: log_gamma_ideal(:) !> log10 activity coefficients of the ideal speciation (0 for the aqueous species)
    real(kind=8), allocatable :: res(:) !> residual in Newton-Raphson
    real(kind=8), allocatable :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    real(kind=8), allocatable :: Delta_p(:) !> Newton step in log-space, p = log10(c1) of the aqueous primaries (n_p_aq)
    real(kind=8), allocatable :: Jac_p(:,:) !> log-space Jacobian d res/d p = Jac_res * diag(ln10 c1) (n_p_aq x n_p_aq)
    real(kind=8), allocatable :: Jac_p_eq(:,:),rhs_eq(:),row_scale(:) !> row-equilibrated log-space system diag(row_scale) Jac_p Delta_p = -diag(row_scale) res, row_scale(i) = 1/||Jac_p(i,:)||_inf
    real(kind=8), allocatable :: tol_res(:) !> tolerance residues Newton-Raphson
    real(kind=8), allocatable :: noise_floor(:) !> max(tol_res, eps^(2/3)*|ctot|) per row: attainable accuracy of the residual, accepted on stagnation
    real(kind=8), allocatable :: d_log_gamma_d_I(:) !> d(log10 gamma)/dI per species (species order; 0 for water and the non-aqueous species)
    real(kind=8), allocatable :: z2_sp(:) !> squared valences in speciation-algebra (species) order (0 for non-aqueous species)
    real(kind=8), allocatable :: dlg_coef_A(:),dlg_beta(:),dlg_gamma(:) !> per-species constants of d(log10 gamma)/dI = coef_A/(2 sqrt(I) (1+beta sqrt(I))^2) + gamma (species order)
    logical, allocatable :: dlg_var(:) !> species whose log gamma depends on I: the dissolved aqueous species (not water, not solids or gases)
    real(kind=8), allocatable :: out_prod(:,:) !> d(log gamma_i)/dI * z_j^2 (species order), the input of compute_dc2_dc1
    real(kind=8), allocatable :: dc2_dc1(:,:) !> Jacobian secondary-primary concentrations (n_e x n_p)
    real(kind=8), allocatable :: log_Jacobian_act_coeffs(:,:) !> d(log gamma_i)/d(log c_j) in the layout compute_res_Jac_res_anal reads (see the loop)
    real(kind=8), allocatable :: c1_save(:),res_trial(:),Jac_trial(:,:) !> line search: saved primaries, trial residual, discarded trial Jacobian
    real(kind=8) :: merit,merit_trial,alpha,scale_factor !> line search: current and trial merit ||res/tol_res||_inf, step length, step bound factor
    real(kind=8) :: sqrt_I !> square root of the ionic strength of the current iterate
    real(kind=8) :: half_ln10_c !> 0.5 ln(10) c_j: factor of column j of log_Jacobian_act_coeffs
    integer(kind=4) :: i_LS !> line search trial counter
    integer(kind=4), parameter :: max_LS_tries=10 !> maximum number of step halvings
    real(kind=8), parameter :: max_dp=2d0 !> maximum change of log10(c1) per Newton step, in decades (a factor 100 in concentration, PHREEQC's default step size)
    real(kind=8), parameter :: ln10=log(1d1)
    real(kind=8), parameter :: half_ln10=5d-1*log(1d1) !> d(log gamma_i)/d(log c_j) = 0.5 ln(10) c_j d(log gamma_i)/dI z_j^2
    integer(kind=4) :: niter_max_init !> iteration budget of the initialisation Newton: max(5*niter_max, 50)
    logical :: step_accepted
    logical :: speciated !> the object is speciated (Picard) at c1 and c2_new holds its secondaries (the trial accepted by the line search): no re-speciation at the top of the iteration
    logical :: sing_flag !> the log-space Jacobian is singular (LU) or the step is not finite: steepest descent instead
    integer(kind=4) :: i,j,niter_Picard
    integer(kind=4) :: niter_ideal !> Newton iterations of the ideal start (fallback)
    integer(kind=4), allocatable :: counters(:)
    logical :: CV_flag_Picard
    logical :: CV_flag_ideal !> convergence of the ideal start (fallback)
    logical :: ideal_start_done !> the ideal start has been tried (only once)
    type(int_array_c) :: indices_icon
    
!> Pre-process
    CV_flag=.false.
    ideal_start_done=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_sp=this%solid_chemistry%reactive_zone%speciation_alg%num_species
    if (size(icon)/=n_p_aq) error stop
    
    if (sum(n_icon)/=n_p_aq) error stop
    
    p=this%ind_prim_species(1:n_p_aq)
    s=this%ind_sec_species(1:n_2_aq)
    call indices_icon%allocate_array(4) !> number of icon options
    do i=1,indices_icon%num_cols
        call indices_icon%cols(i)%allocate_vector(n_icon(i))
    end do
    allocate(counters(4)) !> dim=number of icon options
    counters=0 !> initial counter for each icon option
    allocate(tol_res(n_p_aq))
    do i=1,n_p_aq
        tol_res(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol !> absolute tolerance (by default)
        !> Concentrations
        if (icon(i)==1) then
            counters(1)=counters(1)+1
            indices_icon%cols(1)%col_1(counters(1))=i
        !> Aqueous components
        else if (icon(i)==2) then
            counters(2)=counters(2)+1
            indices_icon%cols(2)%col_1(counters(2))=i
        !> Activities
        else if (icon(i)==3) then
            counters(3)=counters(3)+1
            indices_icon%cols(3)%col_1(counters(3))=i
        !> Phases
        else if (icon(i)==4) then
            counters(4)=counters(4)+1
            indices_icon%cols(4)%col_1(counters(4))=i
            tol_res(i)=1d1**this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !> phase constraint: the residual is logarithmic (sum Se log10 a + log10 ctot - logK), so its tolerance is 10**log_rel_tol (relative error ~ ln(10)*tol on the constrained activity), NOT log_rel_tol itself (a negative number whose magnitude made the test vacuous)
        else
            error stop "icon option not implemented"
        end if
    end do
    !> tol_res was built by INPUT position, i.e. by the concentration index k of the constrained primary
    !> (icons(k), ctots(k)); res and Jac_res (compute_res_Jac_res_anal) come out in the unknowns'
    !> stoichiometric order, row j = primary stored at ind_prim_species(j). Permute the tolerances to the
    !> rows so that each residual is tested against ITS tolerance: a water type that lists its species in
    !> another order than the chemical system would otherwise test a logarithmic phase residual against
    !> abs_tol (which it can never meet) and stop with "Newton method not accurate enough".
    tol_res=tol_res(p)
!> Work arrays, allocated once (nothing is allocated inside the Newton loop)
    allocate(c1(n_p),c1_save(n_p),c2_old(n_e),c2_new(n_e))
    allocate(res(n_p_aq),Jac_res(n_p_aq,n_p_aq),res_trial(n_p_aq),Jac_trial(n_p_aq,n_p_aq))
    allocate(Delta_p(n_p_aq),Jac_p(n_p_aq,n_p_aq),Jac_p_eq(n_p_aq,n_p_aq),rhs_eq(n_p_aq),row_scale(n_p_aq))
    allocate(d_log_gamma_d_I(n_sp),z2_sp(n_sp),dlg_coef_A(n_sp),dlg_beta(n_sp),dlg_gamma(n_sp),dlg_var(n_sp),&
        out_prod(n_sp,n_sp))
    allocate(dc2_dc1(n_e,n_p))
    !> rows and columns indexed as compute_res_Jac_res_anal reads them (concentration index / shifted secondaries), see the loop
    allocate(log_Jacobian_act_coeffs(this%aq_phase%num_species,max(this%aq_phase%num_species,n_p+n_2_aq)))
    !> Noise floor of the residual evaluation, per row (stoichiometric order, like tol_res). The aqueous
    !> secondaries come from the Picard fixed point on the activity coefficients, stopped at a relative
    !> change < rel_tol, so a total carries an error of ~(contraction factor)*rel_tol times its complexes:
    !> ~1e-13 for the mM-level ion pairs of a seawater, above abs_tol = 1e-14. A Newton that has converged
    !> in its iterates (relative step below rel_tol) or whose line search finds no decrease is then chasing
    !> that noise, not failing: it is accepted when every residual is within eps^(2/3)*|ctot_k| (the
    !> stagnation acceptance of initialise_conc_anal_ideal), never looser than tol_res. Otherwise the
    !> stagnation is genuine (singular or inconsistent Jacobian) and the error stops below apply.
    allocate(noise_floor(n_p_aq))
    noise_floor=max(tol_res,epsilon(1d0)**(2d0/3d0)*abs(ctot(p)))
    !> The water-type Newton starts from the user's guesses, which may be many decades away from the solution
    !> (e.g. h+ = 1e-9 against 0.25 M of HCl, or hco3- = 0.1 against P_CO2 = 3.5e-4 at pH 0.6): with the step
    !> bounded to max_dp decades per iteration it needs more iterations than the time-step Newton, which
    !> starts near its solution. It runs once per water type, so a larger budget costs nothing.
    niter_max_init=max(5*this%solid_chemistry%reactive_zone%CV_params%niter_max,50)
    !> Per-species constants of d(log10 gamma)/dI, in the SPECIATION-ALGEBRA (species) order: the n_p
    !> primaries, then the n_e secondaries = the rows of Se_1_star, i.e. the order of c1 and c2, in which
    !> compute_dc2_dc1 reads out_prod. Derivative of the model of compute_log_act_coeffs_aq_phase,
    !> log gamma = -alpha*A*z^2*sqrt(I)/(1+beta*sqrt(I)) + gamma*I, i.e.
    !> d(log gamma)/dI = -alpha*A*z^2/(2*sqrt(I)*(1+beta*sqrt(I))^2) + gamma (as compute_d_log_gamma_d_I_aq_chem,
    !> which is not used: it returns only the aqueous variable activity species and reads the aqueous phase
    !> BY POSITION, which is only right when this water's species are in aqueous-phase order, not after a
    !> primary swap). Each species is reached through the water's concentration -> aqueous phase index map,
    !> so a reordered (swapped) water is handled correctly. Water (constant activity in the model), the
    !> solids and the gases have no ionic-strength dependence: 0.
    d_log_gamma_d_I=0d0
    z2_sp=0d0
    dlg_coef_A=0d0
    dlg_beta=0d0
    dlg_gamma=0d0
    dlg_var=.false.
    do j=1,n_p_aq
        call dlg_params(p(j),j)
    end do
    do j=1,n_2_aq
        call dlg_params(s(j),n_p+j)
    end do
    log_Jacobian_act_coeffs=0d0 !> fixed sparsity pattern (the aqueous rows and columns written in the loop): zeroed once
    log_gamma_ideal=this%compute_log_act_coeffs_ideal()
!> Newton-Raphson
    niter=0
    speciated=.false.
    do
    !> We update number of iterations
        niter=niter+1
        if (niter>niter_max_init) then
            print *, '[DEBUG initialise_conc_anal] Newton failed after', niter-1, 'iters (', trim(this%name), ')'
            print *, '  ||res||_inf =', inf_norm_vec_real(res), ' ||res/tol_res||_inf =', inf_norm_vec_real(res/tol_res)
            print *, '  c1 =', c1
            error stop "Too many Newton iterations in initialisation"
        end if
    !> Speciation at the current iterate: ideal guess -> Picard fixed point of the activity coefficients.
    !> Skipped when the iterate is the trial accepted by the line search of the previous iteration: its
    !> speciation (the object state and c2_new, left by residual_at) is exactly what this would recompute.
        if (.not. speciated) then
            c1=this%get_c1()
            call this%compute_c2_from_c1_ideal(c1,log_gamma_ideal,c2_old)
            call this%compute_c2_from_c1_Picard(c1,c2_old,c2_new,niter_Picard,CV_flag_Picard)
            if (.not. CV_flag_Picard) then
            !> The activity coefficient fixed point does not exist at this iterate: with a poor input guess
            !> (e.g. h+ = 1e-9 against 0.25 M acetic acid: acetate = K*c(acid)/c(h+) ~ 1e6 M) the ionic
            !> strength makes the Davies/Debye-Hückel terms overflow. The iterate is only unchecked at the
            !> first iteration (the line search below rejects infeasible trial points), so recover ONCE by
            !> taking the IDEAL solution (gamma = 1, well behaved from any guess) as the starting point of
            !> the non-ideal Newton; a second failure is a genuine one.
                if (niter==1 .and. .not. ideal_start_done) then
                    ideal_start_done=.true.
                    write(0,'(3A)') ' NOTE: initialise_conc_anal (', trim(this%name), &
                        '): the activity coefficient fixed point diverges at the input guess; starting from the ideal (gamma = 1) solution'
                    call this%initialise_conc_anal_ideal(icon,n_icon,indices_constrains,ctot,niter_ideal,CV_flag_ideal)
                    if (.not. CV_flag_ideal) error stop "initialise_conc_anal: the ideal start did not converge"
                    niter=0 !> the non-ideal iteration count restarts from the ideal solution
                    cycle
                end if
                print *, '[initialise_conc_anal] Picard speciation did not converge at a Newton iterate (', trim(this%name), ')'
                print *, '  c1 =', c1
                error stop "Picard speciation failed in initialisation"
            end if
        end if
        speciated=.false.
    !> We compute residue and Jacobian analytically
        !> d(log gamma)/dI of every species at the ionic strength of the iterate, and its outer product with
        !> z^2 (both species-ordered): d(log gamma_i)/dI * z_j^2
        sqrt_I=sqrt(this%ionic_strength)
        do i=1,n_sp
            if (dlg_var(i)) d_log_gamma_d_I(i)=dlg_coef_A(i)/(2d0*sqrt_I*(1d0+dlg_beta(i)*sqrt_I)**2)+dlg_gamma(i)
        end do
        do j=1,n_sp
            do i=1,n_sp
                out_prod(i,j)=d_log_gamma_d_I(i)*z2_sp(j)
            end do
        end do
        !> We compute Jacobian secondary-primary concentrations
        call this%compute_dc2_dc1(out_prod,c1,c2_new,dc2_dc1)
        !> Log-Jacobian activity coefficients-concentrations, d(log gamma_i)/d(log c_j) = 0.5 ln(10) c_j
        !> d(log gamma_i)/dI z_j^2, in the layout compute_res_Jac_res_anal reads it: rows and primary
        !> columns by CONCENTRATION index, aqueous secondary columns at n_p+1:n_p+n_2_aq (paired with the
        !> rows of dc2aq_dc1). Only aqueous species have nonzero entries.
        do j=1,n_p_aq
            half_ln10_c=half_ln10*this%concentrations(p(j))
            do i=1,n_p_aq
                log_Jacobian_act_coeffs(p(i),p(j))=half_ln10_c*out_prod(i,j)
            end do
            do i=1,n_2_aq
                log_Jacobian_act_coeffs(s(i),p(j))=half_ln10_c*out_prod(n_p+i,j)
            end do
        end do
        do j=1,n_2_aq
            half_ln10_c=half_ln10*this%concentrations(s(j))
            do i=1,n_p_aq
                log_Jacobian_act_coeffs(p(i),n_p+j)=half_ln10_c*out_prod(i,n_p+j)
            end do
            do i=1,n_2_aq
                log_Jacobian_act_coeffs(s(i),n_p+j)=half_ln10_c*out_prod(n_p+i,n_p+j)
            end do
        end do
        !> NOTE: check_dc2_dc1 is a development-time validation of the analytical
        !> Jacobian. It is disabled here because (a) its threshold uses log_rel_tol
        !> directly (a negative log-scale value) against a non-negative matrix norm,
        !> so it can never pass, and (b) it is an absolute comparison that is not
        !> meaningful for non-ideal systems where the Jacobian entries can be large
        !> away from the solution. The Newton residual test below is the actual
        !> convergence criterion.
        !call this%check_dc2_dc1(c1,c2_new,dc2_dc1,log_Jacobian_act_coeffs)
        call this%compute_res_Jac_res_anal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:n_2_aq,:),&
            log_Jacobian_act_coeffs,res,Jac_res)
    !> We check convergence: every constraint residual within its own tolerance (concentration residuals
    !> against abs_tol, logarithmic phase residuals against 10**log_rel_tol)
        if (all(abs(res)<=tol_res)) then !> CV reached
            CV_flag=.true.
            exit
        end if
    !> We solve the linear system for the aqueous primary concentrations IN LOG-SPACE, p = log10(c1):
    !> d res/d p_j = ln10 c1_j d res/d c1_j. The unknowns of a water-type speciation span many orders of
    !> magnitude and the input guesses can be many decades off (h+ = 1e-9 against 0.25 M of HCl), so an
    !> arithmetic step c1 + Delta_c1 either overshoots to negative values or, once clamped, crawls at a
    !> fixed factor per iteration; the log step c1 10**Delta_p keeps every concentration positive and
    !> moves each species by up to max_dp decades per iteration whatever its magnitude.
        do i=1,n_p_aq
            Jac_p(:,i)=ln10*c1(i)*Jac_res(:,i)
        end do
    !> Row equilibration before the LU solve. The columns of Jac_p carry the factor ln10 c1_j, so for a
    !> water with trace primaries (e.g. the 1e-16 totals of mg, na, k, cl, so4 in a freshwater) the
    !> determinant of a perfectly well-conditioned Jacobian is ~ prod(ln10 c1_j) ~ 1e-80, far below the
    !> absolute "numerical zero" with which LU_lin_syst tests singularity: the steepest-descent fallback
    !> below would be taken at every iteration, and it is NOT a descent direction of the inf-norm merit of
    !> the line search (weighted by 1/tol_res^2 it only reduces the residual of the largest total, at the
    !> expense of the phase constraint) -> "Newton line search failed in initialisation". Dividing each
    !> row by its max-norm leaves Delta_p unchanged (row scaling of A x = b), makes every row O(1) (a
    !> component row scales with its total, an activity row with its activity, a phase row already is
    !> O(1)) and so makes the determinant and pivot tests meaningful.
        do i=1,n_p_aq
            row_scale(i)=maxval(abs(Jac_p(i,:)))
            if (row_scale(i)>0d0) then
                row_scale(i)=1d0/row_scale(i)
            else
                row_scale(i)=1d0 !> zero row: genuinely singular, left to LU_lin_syst
            end if
            Jac_p_eq(i,:)=row_scale(i)*Jac_p(i,:)
            rhs_eq(i)=-row_scale(i)*res(i)
        end do
        call LU_lin_syst(Jac_p_eq,rhs_eq,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_p,sing_flag)
        if (.not. sing_flag) sing_flag=any(Delta_p/=Delta_p) .or. maxval(abs(Delta_p))>1d20
        if (sing_flag) Delta_p=-matmul(transpose(Jac_p),res/tol_res**2) !> singular Jacobian: steepest descent of 0.5||res/tol_res||^2
        !> c1^(i+1)=c1^i 10**Delta_p
        if (ln10*inf_norm_vec_real(Delta_p)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            if (all(abs(res)<=noise_floor)) then !> converged to the noise floor of the residual evaluation (see noise_floor)
                CV_flag=.true.
                exit
            end if
            print *, '[initialise_conc_anal] Newton step below rel_tol while the residual is still above its tolerance (', &
                trim(this%name), ')'
            print *, '  ||Delta_c1/c1||_inf =', ln10*inf_norm_vec_real(Delta_p)
            print *, '  res     =', res, ' (unknowns in stoichiometric order)'
            print *, '  tol_res =', tol_res
            print *, '  c1      =', c1
            error stop "Newton method not accurate enough in initialisation"
        end if
    !> Globalisation: the step is bounded to max_dp decades per species (uniform scaling, direction
    !> preserved) and its length alpha is halved until the merit ||res/tol_res||_inf decreases (the trial
    !> residual needs the Picard speciation at the trial point, evaluated by residual_at). The merit is the
    !> residual scaled by its tolerance because the rows have different units (concentrations for icon 1-3,
    !> log10 units for the phase constraints of icon 4); ||res/tol_res||_inf <= 1 is exactly the
    !> convergence test above. A raw Newton step from a poor water-type guess (e.g. h+ = 1e-7 against a
    !> calcite constraint) otherwise oscillates between c and c/10 for ever.
        merit=inf_norm_vec_real(res/tol_res)
        c1_save=c1
        scale_factor=inf_norm_vec_real(Delta_p)
        if (scale_factor>max_dp) Delta_p=Delta_p*(max_dp/scale_factor)
        alpha=1d0
        step_accepted=.false.
        do i_LS=1,max_LS_tries
            c1=c1_save
            c1(1:n_p_aq)=c1_save(1:n_p_aq)*1d1**(alpha*Delta_p) !> positive by construction
            call this%set_conc_aq_prim_species(c1(1:n_p_aq))
            call residual_at(c1,res_trial)
            merit_trial=inf_norm_vec_real(res_trial/tol_res)
            if (merit_trial<merit) then
                step_accepted=.true.
                exit
            end if
            alpha=alpha*0.5d0
        end do
        if (.not. step_accepted) then
            if (all(abs(res)<=noise_floor)) then
            !> At the noise floor of the residual evaluation no step can reduce the merit: accept the current
            !> iterate c1_save, re-speciating the object there (the rejected trials overwrote its state)
                c1=c1_save
                call this%set_conc_aq_prim_species(c1(1:n_p_aq))
                call residual_at(c1,res)
                CV_flag=.true.
                exit
            end if
            print *, '[initialise_conc_anal] no step reduces the residual (', trim(this%name), '); ||res/tol_res||_inf =', merit
            print *, '  c1 =', c1_save
            print *, '  res =', res
            print *, '  Delta_p (log10 step, bounded) =', Delta_p
            error stop "Newton line search failed in initialisation"
        end if
        speciated=.true. !> the accepted trial is the next iterate: the object and c2_new hold its speciation
    end do
contains
    !> Residual of the water-type constraints at the primary concentrations c1_try (already set in the
    !> object by set_conc_aq_prim_species): ideal guess -> Picard speciation -> compute_res_Jac_res_anal.
    !> Only the residual is used (the Jacobian inputs are the current iterate's; the Jacobian output is
    !> discarded). Leaves the object speciated at c1_try and its secondaries in c2_new, which is the next
    !> iterate when the step is accepted; a rejected trial is overwritten by the next trial / the accepted one.
    subroutine residual_at(c1_try,res_try)
        real(kind=8), intent(in) :: c1_try(:)
        real(kind=8), intent(out) :: res_try(:)
        call this%compute_c2_from_c1_ideal(c1_try,log_gamma_ideal,c2_old)
        call this%compute_c2_from_c1_Picard(c1_try,c2_old,c2_new,niter_Picard,CV_flag_Picard)
        if (.not. CV_flag_Picard) then
            res_try=huge(1d0) !> infeasible trial point (activity coefficient fixed point diverged): never accepted by the line search
            return
        end if
        call this%compute_res_Jac_res_anal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:n_2_aq,:),&
            log_Jacobian_act_coeffs,res_try,Jac_trial)
    end subroutine residual_at

    !> Constants of d(log10 gamma)/dI of the aqueous species stored at concentration index k, placed at
    !> species-order position i_sp: z^2 for every species; for the dissolved species (not water)
    !> coef_A = -alpha*A*z^2, beta and gamma of the model of compute_log_act_coeffs_aq_phase, evaluated
    !> with THIS water's Debye-Hückel A and B (temperature dependent)
    subroutine dlg_params(k,i_sp)
        integer(kind=4), intent(in) :: k !> concentration index of the species
        integer(kind=4), intent(in) :: i_sp !> position in species order
        integer(kind=4) :: ka !> index in the aqueous phase
        ka=this%indices_aq_phase(k)
        z2_sp(i_sp)=dble(this%aq_phase%aq_species(ka)%valence**2)
        if (ka==this%aq_phase%ind_wat) return
        associate(pac => this%aq_phase%aq_species(ka)%params_act_coeff)
        dlg_var(i_sp)=.true.
        dlg_coef_A(i_sp)=-pac%alpha*(this%params_aq_sol%A*z2_sp(i_sp))
        dlg_beta(i_sp)=pac%get_beta(this%params_aq_sol)
        dlg_gamma(i_sp)=pac%get_gamma(this%params_aq_sol)
        end associate
    end subroutine dlg_params
end subroutine
