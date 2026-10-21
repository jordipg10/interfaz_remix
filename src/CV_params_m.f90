!> \file CV_params_m.f90
!> \brief Defines convergence parameters for iterative numerical methods.
!> \details This module provides a public type for storing numerical tolerances and control parameters used in iterative algorithms
!>          (e.g., Newton-Raphson, Picard iteration) for chemical speciation and reactive transport simulations.
!>
!> Key convergence criteria (f = component mass-balance residual, u = U·c component totals):
!> - Absolute tolerance: ||f||∞ < abs_tol
!> - Relative tolerance: ||f||∞ < rel_tol·||u|| (Newton solvers) or ||Δc/c||∞ < rel_tol (Picard / initialisation)
!> - Log-absolute tolerance: log₁₀||f||∞ < log_abs_tol
!> - Log-relative tolerance: log₁₀||f||∞ - log₁₀||u|| < log_rel_tol (Newton), or 10**log_rel_tol as the threshold
!>   on the relative change of log₁₀(c) and on logarithmic (mass-action / phase-constraint) residuals
!>
!> Convention: log_abs_tol and log_rel_tol are stored as base-10 exponents (e.g. -14 means 10⁻¹⁴), matching the
!> *_CV_params.dat files; abs_tol and rel_tol are linear values.
!>
!> Physical interpretation:
!> - abs_tol ensures residuals (mass balance errors) are negligible [M/L³]
!> - rel_tol ensures residuals/changes are small relative to the component totals/concentrations [-]
!> - log_abs_tol / log_rel_tol are the same tests in log space (important for trace species)
!> - niter_max prevents infinite loops in non-converging systems
!> - k_div_max limits time step refinement for stability
!>
!> Typical applications:
!> - Newton-Raphson for chemical equilibrium (speciation)
!> - Picard iteration for reactive transport coupling
!> - Time step adaptation for stiff kinetic reactions
module CV_params_m
    implicit none !< Enforce explicit variable declarations for type safety
    !> \struct CV_params_s
    !> \brief Structure for convergence parameters in iterative solvers.
    !> \details Contains absolute, relative, and logarithmic tolerances for convergence checks, along with control parameters
    !>          for iterative algorithms. Used in Newton-Raphson methods for chemical speciation and reactive transport.
    !>
    !> Convergence criteria (f = component residual, u = U·c component totals):
    !>   - Absolute criterion: ||f||∞ < abs_tol [M/L³]
    !>   - Relative criterion: ||f||∞ < rel_tol·max(||u||, ε_mach) (Newton) or ||Δc/c||∞ < rel_tol (Picard) [-]
    !>   - Log-absolute criterion: log₁₀||f||∞ < log_abs_tol (base-10 exponent) [-]
    !>   - Log-relative criterion: log₁₀||f||∞ - log₁₀||u|| < log_rel_tol (Newton), or 10**log_rel_tol as the
    !>     threshold on the relative change of log₁₀(c) and on logarithmic residuals (Picard / initialisation) [-]
    !>
    !> Control parameters:
    !>   - eps: Perturbation of the concentrations for finite-difference Jacobians [M/L³]
    !>   - zero: Threshold below which values treated as zero; also the LU pivot tolerance [M/L³]
    !>   - control_factor: Bound on the per-iteration concentration change ratio, c_new ∈ [f·c_old, c_old/f] (0 < f < 1) [-]
    !>   - niter_max: Safety limit to prevent infinite loops [#]
    !>   - k_div_max: Maximum time step subdivisions for stability [#]
    !>   - est_prm: Time-extrapolation weight of the downstream kinetic reaction rates in the mixing term [-]
    !>
    !> Input: only abs_tol, log_abs_tol, rel_tol, log_rel_tol, est_prm (and optionally niter_max) are read from the
    !> optional <root>_CV_params.dat file (see read_CV_params); eps, zero, control_factor and k_div_max keep their
    !> default initialisation unless set from the driver.
    !>
    !> Physical meaning:
    !> - Smaller tolerances → more accurate but slower convergence
    !> - Logarithmic tolerances important for species spanning many orders of magnitude
    !> - Control factor prevents overshooting in Newton iteration (at most one decade per iteration with 0.1)
    type, public :: CV_params_s !< Public structure holding all convergence parameters for iterative methods
        !> \var abs_tol Arithmetic absolute tolerance for Newton residuals
        !> Convergence criterion: ||f||∞ < abs_tol where f is the component residual (mass balance error)
        !> Default: 1×10⁻¹⁴ [M/L³] - very strict for accurate chemical equilibrium
        !> Physical meaning: maximum acceptable mass balance error in aqueous speciation
        !> Also the absolute local-error tolerance of the RKF45 kinetic integrator
        real(kind=8) :: abs_tol=1d-14 !< [M/L³] absolute tolerance for Newton residuals (mass balance errors)

        !> \var log_abs_tol Logarithmic (base-10 exponent) absolute tolerance for Newton residuals
        !> Convergence criterion: log₁₀||f||∞ < log_abs_tol (the log-space counterpart of abs_tol)
        !> Default: -14 [-], i.e. 10⁻¹⁴ M/L³. Stored as an exponent, matching the *_CV_params.dat files;
        !> read_CV_params warns when it lies below log₁₀(machine epsilon) ≈ -15.65
        real(kind=8) :: log_abs_tol=-14 !< [-] log₁₀ of the absolute tolerance for Newton residuals

        !> \var log_rel_tol Logarithmic (base-10 exponent) relative tolerance
        !> Convergence criteria:
        !>   - Newton solvers: log₁₀||f||∞ - log₁₀ max(||u||, ε_mach) < log_rel_tol (residual relative to the totals)
        !>   - Picard / initialisation: 10**log_rel_tol is the threshold on the relative change of log₁₀(c) and on
        !>     the logarithmic (mass-action, phase-constraint) residuals
        !> Default: -9 [-], i.e. 10⁻⁹ - ensures logarithmic concentration changes negligible
        !> Physical meaning: acceptable change in orders of magnitude for concentrations
        !> Important for trace species where absolute changes are tiny but relative changes significant
        real(kind=8) :: log_rel_tol=-9 !< [-] log₁₀ of the relative tolerance (log-space residuals / concentration changes)
        
        !> \var rel_tol Relative tolerance for convergence
        !> Convergence criteria:
        !>   - Newton solvers: ||f||∞ < rel_tol·max(||u||, ε_mach), the residual relative to the component totals u = U·c
        !>   - Picard / initialisation / surface complexation: ||Δc/c||∞ < rel_tol, the relative solution update
        !>   - RKF45 kinetic integrator: relative local-error tolerance
        !> Default: 1×10⁻⁹ [-]
        !> Physical meaning: fractional residual or fractional change of the solution vector between iterations
        real(kind=8) :: rel_tol=1d-9 !< [-] relative tolerance for residuals / solution updates (dimensionless ratio)
        
        !> \var eps Epsilon for incremental coefficients
        !> Used in finite difference approximations: df/dc ≈ [f(c+ε) - f(c)]/ε
        !> (Jacobians of the incremental-coefficient routines)
        !> Default: 1×10⁻¹² [M/L³] - balance between truncation and roundoff errors
        !> Physical meaning: small additive perturbation of a concentration for numerical Jacobian computation
        !> Too small → roundoff errors dominate; too large → poor derivative approximation
        !> Not read from *_CV_params.dat (default only)
        real(kind=8) :: eps=1d-12 !< [M/L³] perturbation for numerical derivatives (finite differences)
        
        !> \var zero Threshold for zero (may be too high for some applications)
        !> Values |x| < zero treated as zero to avoid division by zero and underflow; also the pivot tolerance
        !> passed to LU_lin_syst, the floor of relative-change denominators and the tolerance of the speciation algebra
        !> Default: 1×10⁻²⁰ [M/L³] - may need reduction for trace species studies
        !> Physical meaning: minimum detectable concentration (practical zero)
        !> Note: May need adjustment for ultra-trace species (e.g., radionuclides)
        !> Not read from *_CV_params.dat (default only)
        real(kind=8) :: zero=1d-20 !< [M/L³] threshold below which values treated as zero
        
        !> \var control_factor Bounds the per-iteration change of the concentrations in the Newton algorithm
        !> Clamp on the update (update_conc_* routines): c_new = min(max(c_old + Δc, f·c_old), c_old/f), i.e. each
        !> concentration may change by at most a factor 1/f per iteration; the full Newton step is applied when it
        !> stays inside these bounds (this is not a step-length multiplier)
        !> Default: 0.1 [-] - at most one order of magnitude up or down per iteration
        !> Range: (0, 1) where smaller values allow larger jumps; f = 1 would freeze the iterate (runtime check in
        !> update_conc_aq_species)
        !> Helps convergence for highly nonlinear systems (e.g., pH near buffering points)
        !> Not read from *_CV_params.dat (default only)
        real(kind=8) :: control_factor=1d-1 !< [-] bound on the concentration change ratio per Newton iteration (0 < f < 1)
        
        !> \var niter_max Maximum number of iterations allowed
        !> Safety limit to prevent infinite loops in non-converging systems (Newton and Picard loops; the
        !> initialisation of the water types uses max(5*niter_max, 50))
        !> Default: 50 [#]. Optional 6th value of the CONVERGENCE PARAMETERS block of *_CV_params.dat
        !> (after est_prm); absent = default. Newton on a stiff kinetic term (rate ∝ 10^p, p = log10 c)
        !> converges linearly far from the root, moving 1/ln10 ≈ 0.43 decades per iteration, so an
        !> injection front several decades off its end-of-step pH needs more than 10 iterations (an
        !> earlier default); a converged loop exits early, so the budget only costs in failing cases.
        !> Physical meaning: maximum attempts before declaring convergence failure
        !> Triggers error or time step reduction if exceeded
        integer(kind=4) :: niter_max=50 !< [#] maximum iterations before convergence failure
        
        !> \var k_div_max Maximum number of time step divisions allowed
        !> Limits time step refinement (halvings of the reactive-mixing sub-step) when convergence fails
        !> Default: 5 [#] - allows Δt to be reduced to Δt/2^5 = Δt/32 for stiff problems
        !> Physical meaning: how many times can halve time step for stiff problems
        !> Prevents excessively small time steps (would make simulation impractically slow); exhausting the
        !> budget is a fatal error ("Too many time step divisions")
        !> Not read from *_CV_params.dat (default only)
        integer(kind=4) :: k_div_max=5 !< [#] maximum time step subdivisions for stability
        
        !> \var est_prm Parameter to estimate kinetic reaction rates of downstream waters
        !> Time extrapolation of the downstream kinetic rates in the reactive mixing term:
        !> r_est = (1 + est_prm)·r_old - est_prm·r_old_old
        !> Default: 0.0 [-] - previous-step rate, no extrapolation (conservative); 1.0 = full linear extrapolation
        !> Physical meaning: weighting factor for estimating the reaction amounts contributed by downstream waters
        !> Fifth value of the CONVERGENCE PARAMETERS block of *_CV_params.dat
        real(kind=8) :: est_prm=0.0 !< [-] time-extrapolation weight of downstream kinetic rates (0 = no extrapolation)
    contains
        !> \brief Reads the convergence parameters from the optional <root>_CV_params.dat file
        procedure :: read_CV_params
    end type CV_params_s

    contains
        !> \brief Reads convergence parameters from the optional <root>_CV_params.dat file
        !> \details The file is optional: when it does not exist the structure is left untouched, so the
        !>          default initialisation of CV_params_s (or whatever the driver set before the call) is used.
        !>          A file that exists but cannot be opened is still a fatal error.
        subroutine read_CV_params(this, path, root)
            class(CV_params_s) :: this
            character(len=*), intent(in) :: path !< Path to input files
            character(len=*), intent(in) :: root !< Root name for input files
            !integer(kind=4), intent(in) :: unit !< Input file unit number
            integer(kind=4) :: unit !< Input file unit number
            integer(kind=4) :: ierr !< I/O error flag
            integer(kind=4) :: ierr_opt !< parse status of the optional niter_max line
            integer(kind=4) :: niter_opt !< niter_max read from the optional 6th line
            logical :: file_exists !< whether <root>_CV_params.dat exists
            character(len=:), allocatable :: fname !< Full path to input file
            character(len=256) :: line !< Record buffer

            fname = trim(path)//trim(root)//'_CV_params.dat'
            inquire(file=fname, exist=file_exists)
            if (.not. file_exists) then
                write(*,*) 'Note: CV_params file not found ('//fname//'): using default convergence parameters'
                return
            end if
            unit = 111 !< Arbitrary unit number for reading CV parameters
            open(unit,file=fname,status='old',action='read',iostat=ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: Could not open CV_params file: '//fname
                error stop
            end if
            do
                read(unit,*,iostat=ierr) line !< Read until end of file
                if (trim(line)=='end') exit
                if (ierr /= 0) exit
                if (trim(line)=='CONVERGENCE PARAMETERS') then
                    read(unit,*) this%abs_tol
                    read(unit,*) this%log_abs_tol
                    read(unit,*) this%rel_tol
                    read(unit,*) this%log_rel_tol
                    read(unit,*) this%est_prm
                    !> Optional 6th value: niter_max (maximum Newton/Picard iterations). The next record is
                    !> either that integer (possibly followed by a comment) or the separator line of the block;
                    !> a record that does not parse as a positive integer leaves the default and is simply skipped
                    read(unit,'(A)',iostat=ierr) line
                    if (ierr==0) then
                        read(line,*,iostat=ierr_opt) niter_opt
                        if (ierr_opt==0 .and. niter_opt>0) this%niter_max=niter_opt
                    end if
                    ! read(unit,*) this%eps
                    ! read(unit,*) this%zero
                    ! read(unit,*) this%control_factor
                    ! read(unit,*) this%k_div_max
                    if (this%abs_tol < epsilon(1d0)) then
                        print *, "WARNING: abs_tol =", this%abs_tol, &
                            " is below machine epsilon =", epsilon(1d0)
                        print *, "  Newton residuals cannot reach this level."
                        print *, "  Consider using abs_tol >= 1d-12."
                    end if
                    if (this%log_abs_tol < log10(epsilon(1d0))) then
                        print *, "WARNING: log_abs_tol =", this%log_abs_tol, &
                            " is below log10(machine epsilon) =", log10(epsilon(1d0))
                        print *, "  Consider using log_abs_tol >= -12."
                    end if
                end if
            end do
            close(unit)
            ! read(unit,*) this%est_prm
        end subroutine read_CV_params
end module CV_params_m