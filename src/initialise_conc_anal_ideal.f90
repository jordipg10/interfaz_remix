!> \file initialise_conc_anal_ideal.f90
!> \brief Newton speciation under ideal conditions (unit activity coefficients).
!> \details
!> Computes species concentrations from water-type constraints using
!> Newton's method with the analytical Jacobian, assuming ideal
!> (unit) activity coefficients.  Used during CHEPROO water-type reading.
!>
!> \author Jordi Petchamé-Guerrero
!> \date Unknown
!> \ingroup chemistry
!> \see initialise_conc_anal, initialise_conc_incr_coeff, aqueous_chemistry_m

!> Computes species concentrations with data from water type definition using Newton method
!> We assume the initial guess of primary aqueous species is already set in the aqueous chemistry object
!> We assume ideal conditions & no exchange reactions
!! We assume concentrations are given in molality
!> This subroutine is only used to read water types based in CHEPROO
subroutine initialise_conc_anal_ideal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag) !< Newton speciation subroutine entry point (ideal activity coefficients)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !< Import aqueous chemistry derived type
    use vectors_m, only: inf_norm_vec_real !< Import infinity-norm utility for real vectors
    use arrays_m, only: int_array_c !< Import integer array container type
    use metodos_sist_lin_m, only: LU_lin_syst !< Import LU-decomposition linear system solver
    implicit none !< Disable implicit typing for safety
!> Arguments
    class(aqueous_chemistry_c) :: this !< Aqueous chemistry object (polymorphic, passed-object dummy)
    integer(kind=4), intent(in) :: icon(:) !> initial condition type
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> indices of constrains in reactive zone stoichiometric matrix
    real(kind=8), intent(in) :: ctot(:) !> data given (must be ordered as in "aq_phase" attribute)
    integer(kind=4), intent(out) :: niter !> number of iterations Newton-Raphson
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    real(kind=8), allocatable :: c1(:),c2(:),log_gamma_vec(:),dc2_dc1(:,:) !< Primary concs, secondary concs, log activity coefficients, Jacobian dc2/dc1
    real(kind=8), allocatable :: res(:) !> residual in Newton-Raphson
    real(kind=8), allocatable :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    real(kind=8), allocatable :: Delta_c1(:),Delta_c1_aux(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: tol_res(:) !> tolerance residues Newton-Raphson
    integer(kind=4) :: i,n_p_aq !< Loop index, number of primary aqueous species
    integer(kind=4), allocatable :: counters(:) !< icon-type counters
    logical :: zero_flag !< Flag for zero concentration detection
    type(int_array_c) :: indices_icon !< Ragged array storing species indices grouped by icon type
    !> Line-search and robustness variables
    real(kind=8) :: res_norm, res_trial_norm, alpha !< Current residual norm, trial residual norm, line-search step length
    real(kind=8), allocatable :: c1_save(:), c2_trial(:) !< Saved primary concs before update, trial secondary concs
    real(kind=8), allocatable :: res_trial(:) !< Trial residual vector after line-search step
    real(kind=8), allocatable :: scaled_Delta(:) !< Step vector scaled by line-search parameter alpha
    logical :: sing_flag !< Flag indicating singular or ill-conditioned Jacobian
    integer(kind=4) :: i_LS !< Line-search backtracking iteration index
    integer(kind=4), parameter :: max_LS = 8 !> max backtracking steps
    !> Levenberg-Marquardt variables
    real(kind=8) :: lambda_LM !< Levenberg-Marquardt damping parameter
    real(kind=8), allocatable :: Jac_LM(:,:) !< Damped Jacobian matrix J + lambda*I
    integer(kind=4) :: i_LM !< Levenberg-Marquardt iteration index
    integer(kind=4), parameter :: max_LM_tries = 10 !< Maximum number of LM damping increases per Newton step
    logical :: step_accepted !< Flag indicating whether LM+line-search produced a descent step
    !> Step bounding variables (concentration-space equivalent of log-space max_dp)
    real(kind=8), parameter :: max_dp = 10d0             !> Max allowed equivalent log-step magnitude per component
    real(kind=8), parameter :: ln10 = log(10d0)          !> Natural logarithm of 10
    real(kind=8) :: scale_factor                          !> Uniform scaling factor for step bounding
    !> Stagnation detection variables
    real(kind=8) :: res_norm_prev                         !> Previous iteration residual norm
    integer(kind=4) :: n_stag                             !> Consecutive stagnation counter
    integer(kind=4), parameter :: n_stag_max = 20         !> Max stagnation iterations before early exit
    real(kind=8), parameter :: stag_rtol = 1d-3           !> Relative tolerance for stagnation detection
    !> Best-solution tracking variables
    real(kind=8), allocatable :: c1_best(:)               !> Best primary concentrations found
    real(kind=8) :: best_res_norm                         !> Smallest residual norm encountered
    !> Perturbation retry variables
    logical :: retry_applied                              !> Whether perturbation retry has been attempted (at most once per call)
    real(kind=8), allocatable :: c1_orig(:)               !> Original initial guess saved for geometric mean retry
    !> Cached loop-invariant parameters (avoid repeated deep struct dereferences per iteration)
    integer(kind=4) :: n_sec_aq                           !< Number of secondary aqueous species
    integer(kind=4) :: niter_max                          !< Maximum Newton iterations
    real(kind=8) :: abs_tol                               !< Absolute convergence tolerance
    real(kind=8) :: rel_tol                               !< Relative convergence tolerance
    real(kind=8) :: log_abs_tol                           !< Log absolute tolerance
    real(kind=8) :: log_rel_tol                           !< Log relative tolerance
    real(kind=8) :: zero_tol                              !< Near-zero threshold
    real(kind=8) :: ctot_norm                             !< ||ctot||_inf (loop-invariant)
    real(kind=8) :: tol_res_norm                          !< ||tol_res||_inf (loop-invariant)
    real(kind=8) :: ctot_scale                            !< max(||ctot||_inf, 1) (loop-invariant)
    real(kind=8) :: eps_d                                 !< Machine epsilon (cached)
    real(kind=8) :: eps23_d                               !< eps^(2/3) — very strict noise floor for initialization
    logical :: skip_eval                                  !< Skip residual/Jacobian recomputation after accepted step
    real(kind=8), allocatable :: U_aq(:,:)                  !< Cached component matrix for aqueous species (loop-invariant)
    integer(kind=4), allocatable :: ind_aq_reorder(:)       !< Cached reordered index array [ind_prim, ind_sec]
    integer(kind=4) :: ind_cstr_LS                          !< Constraint counter for inline residual in LS loop
    
!> Pre-process
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !< Store number of primary aqueous species locally
    CV_flag=.false. !< Initialise convergence flag to false
    if (size(icon)/=n_p_aq) then !< Validate icon array dimension matches number of primary species
        error stop "Dimension error in icon" !< Abort if icon dimension is inconsistent
    else if (sum(n_icon)/=n_p_aq) then !< Validate that n_icon entries sum to total primary species count
         error stop "Dimension error in n_icon" !< Abort if n_icon sum is inconsistent
    end if !< End dimension validation
    !> Cache loop-invariant parameters from deep struct dereferences
    n_sec_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species !< Cache number of secondary aqueous species
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max !< Cache max Newton iterations
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol !< Cache absolute tolerance
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol !< Cache relative tolerance
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol !< Cache log absolute tolerance
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !< Cache log relative tolerance
    zero_tol=this%solid_chemistry%reactive_zone%CV_params%zero !< Cache near-zero threshold
    eps_d=epsilon(1d0) !< Cache machine epsilon
    eps23_d=eps_d**(2d0/3d0) !< Cache eps^(2/3) ≈ 3.65e-11 (very strict noise floor for initialization)
    allocate(c2(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)) !< Allocate secondary species concentration vector
    call indices_icon%allocate_array(4) !> number of icon options
    do i=1,indices_icon%num_cols !< Loop over each icon type to allocate storage
        call indices_icon%cols(i)%allocate_vector(n_icon(i)) !< Allocate vector of size n_icon(i) for icon type i
    end do !< End icon-type allocation loop
    allocate(counters(4)) !> dim=number of icon options
    counters=0 !> initial counter for each icon option
    allocate(tol_res(n_p_aq)) !< Allocate convergence tolerance vector for each primary species equation
    do i=1,n_p_aq !< Loop over primary aqueous species to set tolerances and classify icon types
        tol_res(i)=abs_tol !> absolute tolerance (by default)
        !> Concentrations
        if (icon(i)==1) then !< icon=1: species i is constrained by known concentration
            counters(1)=counters(1)+1 !< Increment concentration-type counter
            indices_icon%cols(1)%col_1(counters(1))=i !< Record species index for concentration constraint
        !> Aqueous components
        else if (icon(i)==2) then !< icon=2: species i is constrained by total component concentration
            counters(2)=counters(2)+1 !< Increment component-type counter
            indices_icon%cols(2)%col_1(counters(2))=i !< Record species index for component constraint
        !> Activities
        else if (icon(i)==3) then !< icon=3: species i is constrained by known activity
            counters(3)=counters(3)+1 !< Increment activity-type counter
            indices_icon%cols(3)%col_1(counters(3))=i !< Record species index for activity constraint
        !> Constrains
        else if (icon(i)==4) then !< icon=4: species i is constrained by mineral equilibrium (saturation index)
            counters(4)=counters(4)+1 !< Increment mineral-constrain counter
            indices_icon%cols(4)%col_1(counters(4))=indices_constrains(counters(4),1) !< Record column index from constrains matrix
            tol_res(i)=abs_tol !> saturation index log10(Omega) is dimensionless, use abs_tol
        else !< Unrecognised icon value
            error stop "icon option not implemented" !< Abort for invalid icon type
        end if !< End icon classification for species i
    end do !< End loop over primary species for icon classification
    !> Cache loop-invariant norms (ctot and tol_res are constant throughout Newton iterations)
    ctot_norm=inf_norm_vec_real(ctot) !< Cache infinity-norm of constraint data vector
    tol_res_norm=inf_norm_vec_real(tol_res) !< Cache infinity-norm of tolerance vector
    ctot_scale=max(ctot_norm,1d0) !< Cache scaled norm for convergence checks
!> Newton-Raphson
    allocate(res(n_p_aq),Jac_res(n_p_aq,n_p_aq)) !< Allocate residual vector and Jacobian matrix
    allocate(Delta_c1(n_p_aq),Delta_c1_aux(n_p_aq)) !< Allocate Newton step vectors (reordered and auxiliary)
    allocate(log_gamma_vec(this%solid_chemistry%reactive_zone%speciation_alg%num_species)) !< Allocate log activity coefficient vector for all species
    allocate(dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,n_p_aq)) !< Allocate Jacobian matrix dc2/dc1
    allocate(c1_save(n_p_aq)) !< Allocate storage for saving c1 before line-search/LM update
    allocate(c2_trial(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)) !< Allocate trial secondary concentrations for line-search
    allocate(res_trial(n_p_aq)) !< Allocate trial residual vector for line-search evaluation
    allocate(scaled_Delta(n_p_aq)) !< Allocate scaled Newton step for line-search
    allocate(Jac_LM(n_p_aq,n_p_aq)) !< Allocate Levenberg-Marquardt damped Jacobian
    allocate(c1_best(n_p_aq)) !< Allocate storage for best primary concentrations found
    allocate(c1_orig(n_p_aq)) !< Allocate storage for original initial guess
    !> Cache U_aq component matrix (loop-invariant — does not change during Newton iterations)
    allocate(U_aq(n_p_aq, this%aq_phase%num_species))
    U_aq(:,1:n_p_aq) = this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq, 1:n_p_aq)
    U_aq(:,n_p_aq+1:this%aq_phase%num_species) = this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq, &
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)
    !> Cache reordered index array for icon-2 residual evaluation
    ind_aq_reorder = [this%ind_prim_species, this%ind_sec_species]
    
    niter=0 !< Initialise Newton iteration counter
    lambda_LM=1d-3 !< Initialise Levenberg-Marquardt damping parameter
    n_stag=0 !< Initialise stagnation counter
    res_norm_prev=huge(1d0) !< Initialise previous residual norm to machine-max (no prior residual)
    best_res_norm=huge(1d0) !< Initialise best residual norm to machine-max
    retry_applied=.false. !< Mark perturbation retry as not yet attempted
    skip_eval=.false. !< No accepted step yet; must compute residual/Jacobian
    call this%set_act_aq_species() !> we set ideal activities of aqueous species
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    log_gamma_vec=this%get_log_gamma() !> we get log_10 activity coefficients of all species
    c1=this%get_c1() !> we get primary aqueous concentrations
    c1_orig=c1 !> save original initial guess for perturbation retry
    do !< Begin main Newton-Raphson iteration loop
        !> We update number of iterations
        niter=niter+1 !< Increment Newton iteration counter
        if (niter>niter_max) then !< Check if maximum iterations exceeded
            if (.not. retry_applied) then !< Check if perturbation retry has not been used yet
                !> Perturbation retry: restart from geometric mean of initial guess and best iterate
                !> c1_retry_j = sqrt(max(c1_orig_j, tiny) * max(c1_best_j, tiny))
                retry_applied = .true. !< Mark perturbation retry as used
                do i = 1, n_p_aq !< Loop over primary species for geometric mean restart
                    c1(i) = sqrt(max(c1_orig(i), tiny(1d0)) * max(c1_best(i), tiny(1d0))) !< Restart c1 as geometric mean of original guess and best iterate
                end do !< End geometric mean restart loop
                call this%set_conc_aq_prim_species(c1) !< Update object with new primary concentrations
                call this%compute_a1() !< Recompute primary activities from restarted c1
                niter = 0 !< Reset iteration counter for fresh retry
                n_stag = 0 !< Reset stagnation counter
                res_norm_prev = huge(1d0) !< Reset previous residual norm
                lambda_LM = 1d-3 !< Reset LM damping parameter
                best_res_norm = huge(1d0) !< Reset best residual norm tracker
                skip_eval = .false. !< Force re-evaluation after retry
                cycle !< Restart Newton loop with perturbed initial guess
            end if !< End retry_applied check
            !> Restore best solution and re-speciate for consistency
            c1=c1_best !< Restore best primary concentrations found so far
            call this%set_conc_aq_prim_species(c1) !< Update object with best primary concentrations
            call this%compute_a1() !< Recompute primary activities for consistency
            call this%compute_c2_from_c1_ideal(c1,log_gamma_vec,c2) !< Recompute secondary concentrations from best c1
            !> Accept if best residual is within conditioning noise: ||r_best|| <= sqrt(eps)*max(||ctot||,1)
            if (best_res_norm<=eps23_d*ctot_scale) then !< Check if best residual is within conditioning noise
                CV_flag=.true. !< Accept convergence: residual is at machine precision
            else !< Best residual exceeds numerical noise threshold
                print *, "[initialise_conc_anal_ideal] Too many Newton iterations, best residual: ", best_res_norm !< Warn user about non-convergence
                CV_flag=.false. !< Mark as not converged
            end if !< End machine-epsilon acceptance check
            exit !< Exit Newton loop after max iterations
        end if !< End max-iteration check block
        if (skip_eval) then !< Reuse residual/Jacobian from the accepted line-search step
            skip_eval = .false. !< Consume the flag; next iteration will recompute unless step is accepted again
        else !< Compute fresh residual and Jacobian
            call this%compute_c2_from_c1_ideal(c1,log_gamma_vec,c2) !< Compute secondary concentrations from current primary using ideal mass-action law
            call this%compute_dc2_dc1_ideal(c1,c2,dc2_dc1) !< Compute analytical Jacobian dc2/dc1 for ideal system
            call this%compute_res_Jac_res_anal_ideal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:&
                n_sec_aq,:),res,Jac_res) !< Assemble residual and Jacobian from icon constraints
        end if !< End skip_eval check
        res_norm = inf_norm_vec_real(res) !< Compute infinity-norm of residual vector
        !> Track best solution: keep the iterate with the smallest residual norm
        if (res_norm<best_res_norm) then !< Check if current residual is new best
            best_res_norm=res_norm !< Update best residual norm
            c1_best=c1 !< Save current primary concentrations as best
        end if !< End best-solution tracking
        !> Check convergence: machine-epsilon safeguard, OR (absolute AND relative) arithmetic, OR logarithmic
        if ((res_norm<=eps_d*ctot_scale) .or. &
            (res_norm<tol_res_norm .and. &
            res_norm<rel_tol*ctot_scale) .or. &
            (res_norm>0d0 .and. &
            (log10(res_norm)<log_abs_tol .and. &
            (log10(res_norm)-log10(ctot_scale))< &
            log_rel_tol))) then !< End multi-criteria convergence test
            CV_flag=.true. !< Mark as converged
            exit !< Exit Newton loop on convergence
        end if !< End convergence check block
        !> Stagnation detection: check if residual is improving sufficiently
        if (res_norm<res_norm_prev*(1d0-stag_rtol)) then !< Check if residual decreased sufficiently (not stagnating)
            n_stag=0 !< Reset stagnation counter on sufficient decrease
            res_norm_prev=res_norm !< Update previous residual norm for next comparison
        else !< Residual did not decrease sufficiently
            n_stag=n_stag+1 !< Increment consecutive stagnation counter
            if (n_stag>=n_stag_max) then !< Check if stagnation limit reached
                if (.not. retry_applied) then !< Check if perturbation retry is still available
                    !> Stagnation detected: try perturbation retry before giving up
                    retry_applied = .true. !< Mark retry as used
                    do i = 1, n_p_aq !< Loop over primary species for geometric mean restart
                        c1(i) = sqrt(max(c1_orig(i), tiny(1d0)) * max(c1_best(i), tiny(1d0))) !< Geometric mean of original and best iterate
                    end do !< End geometric mean restart loop
                    call this%set_conc_aq_prim_species(c1) !< Update object with restarted primary concentrations
                    call this%compute_a1() !< Recompute primary activities
                    niter = 0 !< Reset iteration counter
                    n_stag = 0 !< Reset stagnation counter
                    res_norm_prev = huge(1d0) !< Reset previous residual norm
                    lambda_LM = 1d-3 !< Reset LM damping parameter
                    best_res_norm = huge(1d0) !< Reset best residual norm
                    skip_eval = .false. !< Force re-evaluation after retry
                    cycle !< Restart Newton loop with perturbed initial guess
                end if !< End stagnation retry check
                !> Restore best solution and re-speciate for consistency
                c1=c1_best !< Restore best primary concentrations found
                call this%set_conc_aq_prim_species(c1) !< Update object with best concentrations
                call this%compute_a1() !< Recompute primary activities for consistency
                call this%compute_c2_from_c1_ideal(c1,log_gamma_vec,c2) !< Recompute secondary concentrations from best c1
                !> Accept if best residual is within conditioning noise
                if (best_res_norm<=eps23_d*ctot_scale) then !< Check if residual is within conditioning noise
                    CV_flag=.true. !< Accept convergence
                else !< Best residual exceeds machine-precision noise
                    print *, "[initialise_conc_anal_ideal] Newton stagnated after", niter, "iterations, residual: ", best_res_norm !< Warn about stagnation
                end if !< End stagnation acceptance check
                exit !< Exit Newton loop due to stagnation
            end if !< End n_stag_max check
        end if !< End stagnation else-branch
        !> Save current state for LM + line search
        c1_save = c1 !< Save current primary concentrations before LM/line-search update
        step_accepted = .false. !< Initialise step acceptance flag to false
        !> Levenberg-Marquardt loop: try increasing lambda until step is accepted
        do i_LM = 1, max_LM_tries !< Loop over LM damping levels
            !> Build damped Jacobian: J_LM = J + lambda*I
            Jac_LM = Jac_res !< Copy Jacobian for damping
            do i = 1, n_p_aq !< Loop over diagonal entries
                Jac_LM(i,i) = Jac_LM(i,i) + lambda_LM !< Add LM damping to diagonal
            end do !< End diagonal damping loop
            !> Solve (J + lambda*I) * Delta = -res
            call LU_lin_syst(Jac_LM,-res,abs_tol,Delta_c1_aux,sing_flag) !< Solve damped linear system for Newton step
            if (sing_flag) then !< Check if LU factorisation detected singularity
                lambda_LM = min(lambda_LM * 10d0, 1d8) !< Increase damping to regularise singular system
                cycle !< Retry with larger lambda
            end if !< End singularity check
            !> Check for NaN or huge values in the Newton step
            sing_flag=.false. !< Reset singularity flag for NaN/overflow check
            do i=1,n_p_aq !< Loop over step components to check for invalid values
                if (Delta_c1_aux(i)/=Delta_c1_aux(i) .or. abs(Delta_c1_aux(i))>1d20) then !< Detect NaN (self-inequality) or overflow
                    sing_flag=.true. !< Mark step as invalid
                    exit !< Stop checking remaining components
                end if !< End NaN/overflow check for component i
            end do !< End NaN/overflow validation loop
            if (sing_flag) then !< Check if step contains invalid values
                lambda_LM=min(lambda_LM*10d0,1d8) !< Increase damping to suppress bad step
                cycle !< Retry with larger lambda
            end if !< End invalid-step check
            !> Check if relative update is negligible
            if (inf_norm_vec_real(Delta_c1_aux/max(abs(this%concentrations(1:n_p_aq)),&
                zero_tol)) &
                < rel_tol**2) then !< Relative update below squared tolerance threshold
                CV_flag=.true. !< Mark as converged based on negligible step
                exit !< Exit LM loop on relative-update convergence
            end if !< End relative update check
            !> Reorder Delta to primary species ordering
            Delta_c1=Delta_c1_aux(this%ind_prim_species) !< Reorder auxiliary step to primary species ordering
            !> Bound step: limit max change per component (concentration-space equivalent of log-space max_dp)
            scale_factor=1d0 !< Initialise step scale factor to 1 (no scaling)
            do i=1,n_p_aq !< Loop over primary species to compute bounding scale
                if (abs(Delta_c1(i))>0d0 .and. c1_save(i)>0d0) then !< Only bound non-zero steps for positive concentrations
                    scale_factor=min(scale_factor, max_dp*ln10*c1_save(i)/abs(Delta_c1(i))) !< Tighten scale to keep step within max_dp log-units
                end if !< End per-component bounding check
            end do !< End step bounding loop
            if (scale_factor<1d0) Delta_c1=Delta_c1*scale_factor !< Apply uniform scale factor if step exceeds bound
            !> Backtracking line search (residual-only: Jacobian deferred to after acceptance)
            alpha = 1d0 !< Initialise line-search step length to full Newton step
            do i_LS = 1, max_LS !< Backtracking line-search loop
                c1 = c1_save !< Reset primary concentrations to saved state
                call this%set_conc_aq_prim_species(c1_save) !< Reset object state to saved concentrations
                scaled_Delta = alpha * Delta_c1 !< Scale Newton step by current line-search parameter
                call this%update_conc_aq_prim_species(c1,scaled_Delta,zero_flag) !< Apply scaled step to primary concentrations (with positivity guard)
                call this%compute_a1() !< Recompute primary activities after update
                call this%compute_c2_from_c1_ideal(c1,log_gamma_vec,c2_trial) !< Compute trial secondary concentrations
                !> Inline residual-only evaluation (avoids computing dc2_dc1 and Jac_res per backtrack)
                ind_cstr_LS = 0
                do i = 1, n_icon(1) !< Icon 1: concentration constraints
                    res_trial(indices_icon%cols(1)%col_1(i)) = this%concentrations(indices_icon%cols(1)%col_1(i)) - &
                        ctot(indices_icon%cols(1)%col_1(i))
                end do
                do i = 1, n_icon(2) !< Icon 2: aqueous component constraints
                    res_trial(indices_icon%cols(2)%col_1(i)) = dot_product(U_aq(indices_icon%cols(2)%col_1(i),:), &
                        this%concentrations(ind_aq_reorder)) - ctot(indices_icon%cols(2)%col_1(i))
                end do
                do i = 1, n_icon(3) !< Icon 3: activity constraints
                    res_trial(indices_icon%cols(3)%col_1(i)) = this%activities(indices_icon%cols(3)%col_1(i)) - &
                        ctot(indices_icon%cols(3)%col_1(i))
                end do
                do i = 1, n_icon(4) !< Icon 4: phase equilibrium constraints
                    ind_cstr_LS = ind_cstr_LS + 1
                    res_trial(indices_icon%cols(4)%col_1(i)) = dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se( &
                        indices_constrains(ind_cstr_LS,2), 1:n_p_aq), &
                        log10(this%activities(this%ind_prim_species))) + &
                        log10(ctot(indices_icon%cols(4)%col_1(i))) - &
                        log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr_LS,2))%eq_cst)
                end do
                res_trial_norm = inf_norm_vec_real(res_trial) !< Compute infinity-norm of trial residual
                if (res_trial_norm < res_norm) then !< Check if trial step produces descent (residual decreased)
                    step_accepted = .true. !< Mark step as accepted
                    lambda_LM = max(lambda_LM / 10d0, 1d-12) !< Decrease LM damping on successful step
                    exit !< Exit line-search loop on success
                end if !< End descent check
                alpha = alpha * 0.5d0 !< Halve step length for next backtracking attempt
            end do !< End backtracking line-search loop
            if (step_accepted) exit !< Exit LM loop if line-search found a descent step
            !> Line search failed: restore and increase damping
            c1 = c1_save !< Restore primary concentrations after failed line-search
            call this%set_conc_aq_prim_species(c1_save) !< Reset object state after failed line-search
            lambda_LM = min(lambda_LM * 10d0, 1d8) !< Increase LM damping for next attempt
        end do !< End Levenberg-Marquardt loop
        !> Compute Jacobian once for accepted step (deferred from LS loop) and cache for next iteration
        if (step_accepted) then !< Step was accepted by LM+line-search
            call this%compute_dc2_dc1_ideal(c1,c2_trial,dc2_dc1) !< Compute Jacobian dc2/dc1 at accepted point
            call this%compute_res_Jac_res_anal_ideal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:&
                n_sec_aq,:),res,Jac_res) !< Assemble residual and Jacobian at accepted point
            c2 = c2_trial !< Reuse accepted trial secondary concentrations
            skip_eval = .true. !< Signal next iteration to skip residual/Jacobian computation
        end if !< End accepted-step caching
        if (CV_flag) exit !> relative update convergence detected inside LM loop
        !> If all LM+LS attempts failed, restore state and reset lambda
        if (.not. step_accepted) then !< Check if all LM+line-search attempts failed to find a descent step
            c1 = c1_save !< Restore primary concentrations to pre-LM state
            call this%set_conc_aq_prim_species(c1_save) !< Reset object to pre-LM concentrations
            call this%compute_a1() !< Recompute primary activities for consistency
            lambda_LM = 1d-3 !< Reset LM damping to initial value for next Newton iteration
        end if !< End failed-step restoration
    end do !< End main Newton-Raphson iteration loop
    call this%compute_pH() !< Compute pH from converged hydrogen ion activity
    call this%compute_activities_aq() !< Compute activities of all aqueous species from converged concentrations
    call this%compute_salinity() !> we compute salinity to change units
    call this%compute_alkalinity() !< Compute alkalinity from converged speciation
    deallocate(res,Jac_res,Delta_c1,Delta_c1_aux,log_gamma_vec,dc2_dc1,c1_save,c2_trial,res_trial,scaled_Delta,Jac_LM,c1_best,&
        c1_orig,U_aq,ind_aq_reorder) !< Free all locally allocated working arrays
end subroutine !< End initialise_conc_anal_ideal