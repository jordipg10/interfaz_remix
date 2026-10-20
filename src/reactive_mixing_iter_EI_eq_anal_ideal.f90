!> \file reactive_mixing_iter_EI_eq_anal_ideal.f90
!> \brief Reactive mixing iteration with Euler-implicit equilibrium reactions for ideal systems.
!> \details This subroutine performs a single reactive mixing iteration of the WMA (Water Mixing Approach)
!> for equilibrium-only systems (no kinetic reactions). It solves for the variable activity species
!> concentrations \f$ \mathbf{c}_{v} \f$ and component concentrations \f$ \mathbf{u} \f$ after reactive mixing.
!>
!> The approach uses:
!> - **Analytical Jacobians** for the Newton method
!> - **Ideal solution** assumption (unit activity coefficients)
!> - **Equilibrium-only** speciation (no kinetic source terms)
!>
!> **Retry strategy:**
!>
!> The equilibrium problem \f$ \mathbf{U}\,\mathbf{c}_v(\mathbf{c}_1) = \hat{\mathbf{u}} \f$ has no time dependence
!> (no kinetic source term), so sub-stepping in \f$ \Delta t \f$ cannot change it; the only retry lever is
!> the initial guess. If Newton fails, the following guesses are tried in order, each at most once:
!> 1. closest downstream target water (if `c1_downstream` is given),
!> 2. linear extrapolation with \f$ \mu \in \{0, 0.25, 0.5, 0.75, 1\} \f$,
!> 3. the mixed concentrations \f$ \hat{\mathbf{c}} \f$.
!> If all of them fail the run stops with diagnostics (fatal error).
!>
!> **Initial guess strategy:**
!>
!> The extrapolated initial guess for primary species is
!> \f[
!>   c_{1,j}^{(0)} = (1 + \mu)\,c_{1,j}^{(k)} - \mu\,c_{1,j}^{(k-1)}
!> \f]
!> clamped to \f$ c_{1,j}^{(k)} \f$ if the extrapolation would go negative. A concentration floor is
!> applied to every guess: if \f$ c_{1,j}^{(0)} < \max(10^{-20},\,10^{-8}\hat{c}_j) \f$ while
!> \f$ \hat{c}_j > 10^{-20} \f$, the mixed value \f$ \hat{c}_j \f$ is used instead (`rel_floor`; it only
!> replaces guesses that are negligible against the mixed value, not legitimate low concentrations of a
!> parcel that mixes with much richer waters).
!>
!> \param[in,out] this Aqueous chemistry object at current time step; on exit its variable activity
!>                     concentrations, aqueous activities and pH are those of the solution
!> \param[in] c1_old Primary species concentrations at previous time step (n_p)
!> \param[in] c_hat Variable activity species concentrations after mixing (n_v)
!> \param[in] mix_ratio_r Reaction mixing ratio of this target water \f$ \lambda_{r} \f$ [-] (unused: equilibrium only)
!> \param[in] Delta_t Time step \f$ \Delta t \f$ [T] (unused: the equilibrium problem is time independent)
!> \param[in] theta Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-] (unused: equilibrium only)
!> \param[in,out] conc_nc Variable activity species concentrations (already allocated, size n_v); current
!>                        concentrations \f$ \mathbf{c}_v^{(k)} \f$ on entry, solution on exit
!> \param[in,out] conc_comp Component concentrations (already allocated, size n_p); overwritten on exit with \f$ \mathbf{U}\mathbf{c}_{v} \f$
!> \param[in] c1_downstream Optional primary concentrations of the closest downstream target water (n_p), used as first initial guess
!>
!> \pre The aqueous chemistry object must have valid speciation algebra, solid chemistry, and convergence parameters.
!> \pre `conc_nc` must be pre-allocated with size \f$ n_v \f$ and contain physically meaningful values.
!> \post On exit, `conc_nc` contains the converged species concentrations and `conc_comp` contains the
!>       corresponding component concentrations.
!>
!> \note The chemistry object is not updated between Newton attempts: `Newton_EI_eq_ideal` reads only
!>       `u_hat` and `conc_nc` (plus the object's equilibrium constants), so a single update after
!>       convergence suffices.
!>
!> \sa Newton_EI_eq_ideal, initialise_iterative_method

subroutine reactive_mixing_iter_EI_eq_anal_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp,c1_downstream)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method  !< Import aqueous chemistry class and initial guess helper
    implicit none                                                            !< Require explicit declaration of all variables
!> Arguments
    class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: c1_old(:)              !< Primary species concentrations at previous time step \f$ \mathbf{c}_1^{(k-1)} \f$ (n_p)
    real(kind=8), intent(in) :: c_hat(:)               !< Variable activity species concentrations after mixing \f$ \hat{\mathbf{c}}_v \f$ (n_v)
    real(kind=8), intent(in) :: mix_ratio_r            !< Reaction mixing ratio of this target water \f$ \lambda_{r} \f$ [-] (unused: equilibrium only)
    real(kind=8), intent(in) :: Delta_t                !< Time step \f$ \Delta t \f$ [T] (unused: equilibrium problem is time independent)
    real(kind=8), intent(in) :: theta                  !< Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-] (unused: equilibrium only)
    real(kind=8), intent(inout) :: conc_nc(:)          !< Variable activity species concentrations (already allocated) (n_v)
    real(kind=8), intent(inout) :: conc_comp(:)        !< Component concentrations (already allocated) (n_p)
    real(kind=8), intent(in), optional :: c1_downstream(:) !< Primary concentrations of closest downstream target water (n_p)
!> Local variables
    real(kind=8), allocatable :: c1(:)                 !< Primary species concentrations at current time step \f$ \mathbf{c}_1^{(k)} \f$ (n_p)
    real(kind=8), allocatable :: u_hat(:)              !< Component concentrations after mixing \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$ (n_p)
    integer(kind=4) :: n_p                             !< Number of primary species [-]
    integer(kind=4) :: niter                           !< Number of Newton iterations in last Newton call [-]
    integer(kind=4) :: j_sp                            !< Loop index over primary species for concentration floor [-]
    real(kind=8) :: mu                                 !< Newton initialisation parameter \f$ \mu \in [0,1] \f$ [-]
    real(kind=8), parameter :: conc_floor = 1d-20      !< Concentration floor: threshold below which a species is considered effectively zero [-]
    real(kind=8), parameter :: rel_floor = 1d-8        !< Relative floor of the initial guess: a guess below rel_floor * c_hat is replaced by the mixed value [-]
    logical :: CV_flag                                 !< Convergence flag from Newton solver: .true. if converged [-]
    logical :: use_downstream                          !< Next attempt uses the downstream-water initial guess [-]
    logical :: c_hat_tried                             !< Whether the c_hat-based initial guess has been attempted as a last resort [-]

    associate(comp_mat => this%solid_chemistry%reactive_zone%speciation_alg%comp_mat)  !< \f$ \mathbf{U} \f$ (n_p x n_v) by reference: no per-call copy
!> -----------------------------------------------------------------------
!> \section rmix_eq_preproc Pre-processing: dimensions, mixed components, snapshot of current c1
!> -----------------------------------------------------------------------
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species    !< Extract number of primary species from speciation algebra
    allocate(u_hat(n_p),c1(n_p))                                             !< Allocate mixed component vector and primary concentration snapshot
    c1=conc_nc(1:n_p)                                                        !< Snapshot \f$ \mathbf{c}_1^{(k)} \f$ before conc_nc is overwritten by the guesses
    u_hat=matmul(comp_mat,c_hat)                                             !< \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$ (independent of Delta_t: no kinetic term)

!> -----------------------------------------------------------------------
!> \section rmix_eq_process Process: Newton solve over a sequence of distinct initial guesses
!> -----------------------------------------------------------------------
    !> Order of the guesses: downstream water (if given), extrapolation with mu = 0, 0.25, ..., 1,
    !> mixed concentrations c_hat. Each is tried once; sub-stepping in Delta_t would only repeat them.
    use_downstream=present(c1_downstream)                                    !< Downstream guess goes first when available
    mu=-0.25d0                                                               !< First extrapolation attempt uses mu = 0
    c_hat_tried=.false.                                                      !< c_hat fallback not yet attempted
    try_loop: do                                                             !< One Newton attempt per initial guess
        if (use_downstream) then                                             !< Guess 1: closest downstream target water
            use_downstream=.false.                                           !< Consume the option
            conc_nc(1:n_p)=c1_downstream(1:n_p)                              !< Use downstream target water concentrations as initial guess
        else if (mu<1d0) then                                                !< Guess 2: extrapolation, mu = 0, 0.25, ..., 1 (exact in binary)
            mu=mu+0.25d0                                                     !< Next mu
            call initialise_iterative_method(c1_old,c1,mu,conc_nc(1:n_p))   !< \f$ c_{1,j}^{(0)} = (1+\mu)\,c_{1,j}^{(k)} - \mu\,c_{1,j}^{(k-1)} \f$
        else if (.not. c_hat_tried) then                                     !< Guess 3: mixed concentrations
            c_hat_tried=.true.                                               !< Consume the option
            conc_nc(1:n_p)=c_hat(1:n_p)                                      !< Use mixed concentrations as initial guess
        else                                                                 !< Every initial guess has failed
            exit try_loop                                                    !< Leave with CV_flag = .false.
        end if                                                               !< End initial guess selection
        !> Floor near-zero species: replace with mixed values when the guess is negligible
        !> (absolute floor conc_floor and relative floor rel_floor * c_hat)
        do j_sp=1,n_p                                                        !< Loop over each primary species
            if (conc_nc(j_sp) < max(conc_floor, c_hat(j_sp)*rel_floor) .and. & !< Initial guess below absolute or relative floor?
                c_hat(j_sp) > conc_floor) then                                !< Mixed value is not negligible?
                conc_nc(j_sp)=c_hat(j_sp)                                     !< Use mixed concentration as initial guess
            end if                                                            !< End floor check
        end do                                                                !< End initial guess floor loop
        call this%Newton_EI_eq_ideal(u_hat,conc_nc,niter,CV_flag)            !< Reads only u_hat and conc_nc: the object needs no update before the call
        if (CV_flag) exit try_loop                                           !< Converged: leave the retry loop
    end do try_loop                                                          !< End initial guess loop
    if (.not. CV_flag) then                                                  !< Unrecoverable non-convergence
        !> conc_nc holds the best iterate of the last Newton attempt
        print *, "FATAL: reactive_mixing_iter_EI_eq_anal_ideal: Newton failed for every initial guess"
        print *, "  best ||U c - u_hat||_inf=", maxval(abs(matmul(comp_mat,conc_nc)-u_hat)), &
            " ||U c||_inf=", maxval(abs(matmul(comp_mat,conc_nc))), &
            " abs_tol=", this%solid_chemistry%reactive_zone%CV_params%abs_tol, &
            " rel_tol=", this%solid_chemistry%reactive_zone%CV_params%rel_tol
        print *, "  c1_old=", c1_old                                          !< Print old primary concentrations
        print *, "  c_hat=", c_hat                                            !< Print mixed species concentrations
        print *, "  conc_nc=", conc_nc(1:n_p)                                 !< Print best primary concentrations
        error stop "Newton did not converge for any initial guess in reactive_mixing_iter_EI_eq_anal_ideal. Check CV parameters."
    end if                                                                    !< End failure handling

!> -----------------------------------------------------------------------
!> \section rmix_eq_postproc Post-processing: update object, compute component concentrations
!> -----------------------------------------------------------------------
    call this%set_conc_var_act_species(conc_nc)                              !< Set variable activity species concentrations in the chemistry object
    call this%set_act_aq_species()                                            !< Update aqueous species activities based on converged concentrations
    call this%compute_pH()                                                   !< Compute pH from the converged aqueous species concentrations
    conc_comp=matmul(comp_mat,conc_nc)                                       !< Compute final component concentrations \f$ \mathbf{u} = \mathbf{U}\mathbf{c}_v \f$
    end associate                                                            !< End comp_mat association
    deallocate(u_hat,c1)                                                     !< Free locally allocated arrays
end subroutine                                                               !< End reactive_mixing_iter_EI_eq_anal_ideal