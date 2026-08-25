!> \file reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2.f90
!> \brief Reactive mixing iteration with Euler-implicit kinetic reactions for ideal equilibrium-kinetic systems (option 2)
!> \details This subroutine performs a single reactive mixing iteration of the WMA (Water Mixing Approach)
!> using Euler-implicit time integration for kinetic reactions. It solves for the variable activity species
!> concentrations \f$ \mathbf{c}_{v} \f$ and component concentrations \f$ \mathbf{u} \f$ after reactive mixing.
!>
!> The approach uses:
!> - **Analytical Jacobians** for the Newton method
!> - **Ideal solution** assumption (unit activity coefficients)
!> - **Option 2**: averaging kinetic reaction rates between old and new time steps for improved stability
!>
!> **Adaptive time-stepping strategy:**
!>
!> If Newton does not converge, the subroutine first sweeps the initialisation parameter
!> \f$ \mu \in \{0, 0.25, 0.5, 0.75, 1.0\} \f$ to vary the initial guess. If all \f$ \mu \f$ values
!> fail, the time step is halved (\f$ \Delta t \leftarrow \Delta t / 2 \f$) and the process restarts.
!> This continues until \f$ 2^{k_{\text{div}}} \f$ consecutive sub-steps converge, covering the
!> original \f$ \Delta t \f$, or until the maximum number of divisions \f$ k_{\text{div,max}} \f$
!> is exceeded (fatal error).
!>
!> **Initial guess strategy:**
!>
!> The initial guess for primary species is computed via linear extrapolation:
!> \f[
!>   c_{1,j}^{(0)} = (1 + \mu)\,c_{1,j}^{(k)} - \mu\,c_{1,j}^{(k-1)}
!> \f]
!> clamped to \f$ c_{1,j}^{(k)} \f$ if the extrapolation would go negative. A concentration floor
!> \f$ 10^{-20} \f$ is applied: if the initial guess is below this floor but the mixed concentration
!> \f$ \hat{c}_j \f$ is above it, the mixed value is used instead.
!>
!> \param[in,out] this Aqueous chemistry object at current time step
!> \param[in] c1_old Primary species concentrations at previous time step (n_p)
!> \param[in] c_hat Variable activity species concentrations after mixing (n_v)
!> \param[in] mix_ratio_r_old Mixing ratio of kinetic reaction amounts from previous time step \f$ \lambda_{r}^{\mathrm{old}} \f$ [-]
!> \param[in] mix_ratio_r_new Mixing ratio of kinetic reaction amounts at current time step \f$ \lambda_{r}^{\mathrm{new}} \f$ [-]
!> \param[in] Delta_t Time step \f$ \Delta t \f$ [T]
!> \param[in] theta Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-] (0 = explicit, 1 = fully implicit)
!> \param[in,out] conc_nc Variable activity species concentrations (already allocated, size n_v); initial guess on entry, solution on exit
!> \param[in,out] conc_comp Component concentrations (already allocated, size n_p); overwritten on exit with \f$ \mathbf{U}\mathbf{c}_{v} \f$
!>
!> \pre The aqueous chemistry object must have valid speciation algebra, solid chemistry, and convergence parameters.
!> \pre `conc_nc` must be pre-allocated with size \f$ n_v \f$ and contain physically meaningful values.
!> \post On exit, `conc_nc` contains the converged species concentrations and `conc_comp` contains the
!>       corresponding component concentrations.
!>
!> \sa Newton_EI_eq_kin_anal_ideal_opt2, initialise_iterative_method

subroutine reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2(this,c1_old,c_hat,mix_ratio_r,&
    Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method  !< Import aqueous chemistry class and initial guess helper
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: c1_old(:)              !< Primary species concentrations at previous time step \f$ \mathbf{c}_1^{(k-1)} \f$ (n_p)
    real(kind=8), intent(in) :: c_hat(:)               !< Variable activity species concentrations after mixing \f$ \hat{\mathbf{c}}_v \f$ (n_v)
    real(kind=8), intent(in) :: mix_ratio_r        !< Mixing ratio of reaction amounts \f$ \lambda_{r}^{\mathrm{new}} \f$ [-]
    real(kind=8), intent(in) :: Delta_t                !< Time step \f$ \Delta t \f$ [T]
    real(kind=8), intent(in) :: theta                  !< Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-]
    real(kind=8), intent(inout) :: conc_nc(:)          !< Variable activity species concentrations (already allocated) (n_v)
    real(kind=8), intent(inout) :: conc_comp(:)        !< Component concentrations (already allocated) (n_p)
!> Local variables
    real(kind=8), allocatable :: c1(:)                 !< Primary species concentrations at current time step \f$ \mathbf{c}_1^{(k)} \f$ (n_p)
    real(kind=8), allocatable :: u_hat(:)              !< Component concentrations after mixing \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$ (n_p)
    real(kind=8), allocatable :: comp_mat(:,:)         !< Cached component matrix \f$ \mathbf{U} \f$ (n_p x n_v)
    integer(kind=4) :: n_p                             !< Number of primary species [-]
    integer(kind=4) :: k_div_max                       !< Cached maximum allowed time step divisions [-]
    integer(kind=4) :: k_div                           !< Counter of time step halvings performed [-]
    integer(kind=4) :: k                               !< Counter of completed converged sub-steps [-]
    integer(kind=4) :: niter                           !< Number of Newton iterations in last Newton call [-]
    integer(kind=4) :: j_sp                            !< Loop index over primary species for concentration floor [-]
    real(kind=8) :: mu                                 !< Newton initialisation parameter \f$ \mu \in [0,1] \f$ [-]
    real(kind=8) :: Delta_t_bis                        !< Current (possibly halved) sub-step size \f$ \Delta t / 2^{k_{\text{div}}} \f$ [T]
    real(kind=8), parameter :: conc_floor = 1d-20      !< Concentration floor: threshold below which a species is considered effectively zero [-]
    logical :: CV_flag                                 !< Convergence flag from Newton solver: .true. if converged [-]
    logical :: c_hat_fallback_tried                     !< Whether c_hat-based initial guess has been attempted as a last resort [-]

!> -----------------------------------------------------------------------
!> \section rmix_preproc Pre-processing: extract dimensions, cache matrices, compute initial guess
!> -----------------------------------------------------------------------
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species   !< Extract number of primary species from speciation algebra
    comp_mat=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat      !< Cache component matrix \f$ \mathbf{U} \f$ (avoids repeated deep dereference)
    k_div_max=this%solid_chemistry%reactive_zone%CV_params%k_div_max         !< Cache maximum allowed time step divisions
    allocate(u_hat(n_p),c1(n_p))                                             !< Allocate component vector and primary concentration snapshot
    c1=conc_nc(1:n_p)                                                        !< Snapshot current primary concentrations \f$ \mathbf{c}_1^{(k)} \f$
    mu=0d0                                                                   !< Initialise extrapolation parameter: \f$ \mu = 0 \f$ (no extrapolation)
    u_hat=matmul(comp_mat,c_hat)                                             !< Compute mixed component concentrations \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$

!> -----------------------------------------------------------------------
!> \section rmix_process Process: Newton iteration with adaptive time-stepping
!> -----------------------------------------------------------------------
    !> Compute initial guess for primary species by linear extrapolation
    call initialise_iterative_method(c1_old,c1,mu,conc_nc(1:n_p))           !< \f$ c_{1,j}^{(0)} = (1+\mu)\,c_{1,j}^{(k)} - \mu\,c_{1,j}^{(k-1)} \f$
    !> Floor near-zero species: replace with mixed values when initial guess is negligible
    !> Uses both an absolute floor (conc_floor) and a relative floor (1e-4 * c_hat)
    !> to catch cases where the initial guess is orders of magnitude below c_hat
    !> (e.g. sharp concentration fronts arriving at clean water)
    do j_sp=1,n_p                                                            !< Loop over each primary species
        if (conc_nc(j_sp) < max(conc_floor, c_hat(j_sp)*1d-4) .and. &        !< Initial guess below absolute or relative floor?
            c_hat(j_sp) > conc_floor) then                                    !< Mixed value is not negligible?
            conc_nc(j_sp)=c_hat(j_sp)                                         !< Use mixed concentration as initial guess
        end if                                                                !< End floor check
    end do                                                                    !< End initial guess floor loop
    k=0                                                                      !< Reset completed sub-step counter
    k_div=0                                                                  !< Reset time step division counter
    c_hat_fallback_tried=.false.                                              !< c_hat fallback not yet attempted
    Delta_t_bis=Delta_t                                                      !< Start with the full time step
    do                                                                       !< Outer loop: accumulate \f$ 2^{k_{\text{div}}} \f$ converged sub-steps
        do                                                                   !< Inner loop: retry Newton with varying \f$ \mu \f$ or halved \f$ \Delta t \f$
            !> Solve the nonlinear reactive system via Newton-Raphson
            call this%Newton_EI_eq_kin_anal_ideal_opt2(u_hat,mix_ratio_r,Delta_t_bis,theta,conc_nc,niter,& !< Call Newton solver
                CV_flag)                                                      !< Returns convergence flag
            if (CV_flag) then                                                !< Did Newton converge?
                k=k+1                                                         !< Increment completed sub-step counter
                exit                                                          !< Exit inner retry loop — proceed to outer loop check
            end if                                                            !< End convergence check
            !> Newton did not converge: adjust initial guess or halve time step
            if (mu<1d0) then                                                 !< Have we exhausted all \f$ \mu \f$ values?
                mu=mu+0.25d0                                                  !< Increment \f$ \mu \f$: try a more extrapolated initial guess
            else                                                              !< All \f$ \mu \f$ values failed
                mu=0d0                                                        !< Reset \f$ \mu \f$ for next round of attempts
                k_div=k_div+1                                                 !< Increment time step division counter
                if (k_div>k_div_max) then                                     !< Exceeded maximum allowed divisions?
                    if (.not. c_hat_fallback_tried) then                       !< c_hat fallback not yet attempted?
                        !> Last-resort fallback: use c_hat(1:n_p) as initial guess.
                        !> This handles sharp concentration fronts where extrapolation-based
                        !> guesses are orders of magnitude from the solution.
                        c_hat_fallback_tried = .true.                          !< Mark fallback as attempted
                        k_div = 0                                              !< Reset time step division counter
                        k = 0                                                  !< Reset completed sub-step counter
                        Delta_t_bis = Delta_t                                  !< Restore original time step
                        mu = 0d0                                               !< Reset extrapolation parameter
                        conc_nc(1:n_p) = c_hat(1:n_p)                          !< Use mixed concentrations as initial guess
                        cycle                                                  !< Retry inner loop with c_hat-based guess
                    end if                                                     !< End c_hat fallback check
                    print *, "FATAL: k_div_max exceeded. k_div=", k_div, " Delta_t_bis=", Delta_t_bis  !< Report failure diagnostics
                    print *, "  c1_old=", c1_old                              !< Print old primary concentrations
                    print *, "  c_hat=", c_hat                                !< Print mixed species concentrations
                    print *, "  conc_nc=", conc_nc(1:n_p)                     !< Print current primary concentration guess
                    print *, "  rk_old=", this%get_rk_old()                   !< Print old kinetic reaction rates
                    error stop "Too many time step divisions in reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2. &
                        Check CV parameters."                                 !< Fatal error: unrecoverable non-convergence
                end if                                                        !< End k_div_max check
                Delta_t_bis=Delta_t_bis/2d0                                   !< Halve the sub-step size \f$ \Delta t_{\text{bis}} \leftarrow \Delta t_{\text{bis}} / 2 \f$
            end if                                                            !< End mu/k_div branching
            !> Recompute initial guess with updated \f$ \mu \f$ for the next Newton attempt
            call initialise_iterative_method(c1_old,c1,mu,conc_nc(1:n_p))    !< Linear extrapolation with new \f$ \mu \f$
            !> Re-apply concentration floor (absolute and relative)
            do j_sp=1,n_p                                                     !< Loop over each primary species
                if (conc_nc(j_sp) < max(conc_floor, c_hat(j_sp)*1d-4) .and. & !< Below absolute or relative floor?
                    c_hat(j_sp) > conc_floor) then                             !< Mixed value is not negligible?
                    conc_nc(j_sp)=c_hat(j_sp)                                 !< Replace with mixed concentration
                end if                                                        !< End floor check
            end do                                                            !< End retry floor loop
        end do                                                                !< End inner retry loop
        !> One converged solve at the (possibly subdivided) Delta_t_bis IS the result: the former
        !> 2^k_div re-solves repeated the identical (u_hat,Delta_t_bis) system from the converged
        !> state (conc_nc unchanged) and only caused the exponential-cost hang.
        exit                                                                 !< Converged sub-step obtained; exit outer loop
    end do                                                                    !< End outer sub-step accumulation loop

!> -----------------------------------------------------------------------
!> \section rmix_postproc Post-processing: compute component concentrations and deallocate
!> -----------------------------------------------------------------------
    conc_comp=matmul(comp_mat,conc_nc)                                       !< Compute final component concentrations \f$ \mathbf{u} = \mathbf{U}\mathbf{c}_v \f$
    deallocate(u_hat,c1,comp_mat)                                            !< Free all locally allocated arrays
end subroutine                                                               !< End reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2