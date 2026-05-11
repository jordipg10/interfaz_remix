!> \file CV_params_m.f90
!> \brief Defines convergence parameters for iterative numerical methods.
!> \details This module provides a public type for storing numerical tolerances and control parameters used in iterative algorithms
!>          (e.g., Newton-Raphson, Picard iteration) for chemical speciation and reactive transport simulations.
!>
!> Key convergence criteria:
!> - Absolute tolerance: |f(x)| < abs_tol for residual functions
!> - Relative tolerance: |Δx|/|x| < rel_tol for solution updates
!> - Logarithmic tolerance: |log₁₀(c) - log₁₀(c_old)| < log_rel_tol for concentrations
!>
!> Physical interpretation:
!> - abs_tol ensures residuals (mass balance errors) are negligible [M/L³]
!> - rel_tol ensures relative changes in solution are small [-]
!> - log_rel_tol handles logarithmic concentration changes (important for trace species)
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
    !> Convergence criteria:
    !>   - Absolute criterion: ||f(x)|| < abs_tol where f is residual function [M/L³]
    !>   - Relative criterion: ||Δx||/||x|| < rel_tol where Δx is update [-]
    !>   - Logarithmic criterion: |log₁₀(c_new) - log₁₀(c_old)| < log_rel_tol for concentrations [-]
    !>
    !> Control parameters:
    !>   - eps: Small number for numerical derivatives (finite differences) [-]
    !>   - zero: Threshold below which values treated as zero [M/L³]
    !>   - control_factor: Damping factor for Newton updates (0 < factor ≤ 1) [-]
    !>   - niter_max: Safety limit to prevent infinite loops [#]
    !>   - k_div_max: Maximum time step subdivisions for stability [#]
    !>   - est_prm: Estimation parameter for kinetic reaction amounts downstream [-]
    !>
    !> Physical meaning:
    !> - Smaller tolerances → more accurate but slower convergence
    !> - Logarithmic tolerance important for species spanning many orders of magnitude
    !> - Control factor prevents oscillations in Newton iteration (damping)
    type, public :: CV_params_s !< Public structure holding all convergence parameters for iterative methods
        !> \var abs_tol Arithmetic absolute tolerance for Newton residuals
        !> Convergence criterion: ||f(x)|| < abs_tol where f is residual (mass balance error)
        !> Default: 1×10⁻¹⁴ [M/L³] - very strict for accurate chemical equilibrium
        !> Physical meaning: maximum acceptable mass balance error in aqueous speciation
        real(kind=8) :: abs_tol !< [M/L³] absolute tolerance for Newton residuals (mass balance errors)
        real(kind=8) :: log_abs_tol !< [M/L³] absolute tolerance for Picard residuals (mass balance errors)
        !> \var log_rel_tol Logarithmic relative tolerance for Newton residuals
        !> Convergence criterion: |log₁₀(Δc)/log₁₀(c)| < log_rel_tol for concentrations
        !> Default: 1×10⁻⁹ [-] - ensures logarithmic concentration changes negligible
        !> Physical meaning: acceptable change in orders of magnitude for concentrations
        !> Important for trace species where absolute changes are tiny but relative changes significant
        real(kind=8) :: log_rel_tol !< [-] logarithmic relative tolerance for concentration changes
        
        !> \var rel_tol Relative tolerance for convergence
        !> Convergence criterion: ||Δx||/||x|| < rel_tol where Δx is Newton update
        !> Default: 1×10⁻¹⁶ [-] - extremely strict relative convergence
        !> Physical meaning: fractional change in solution vector between iterations
        real(kind=8) :: rel_tol !< [-] relative tolerance for solution updates (dimensionless ratio)
        
        !> \var eps Epsilon for incremental coefficients
        !> Used in finite difference approximations: df/dx ≈ [f(x+ε) - f(x)]/ε
        !> Default: 1×10⁻¹² [-] - balance between truncation and roundoff errors
        !> Physical meaning: small perturbation for numerical Jacobian computation
        !> Too small → roundoff errors dominate; too large → poor derivative approximation
        real(kind=8) :: eps=1d-12 !< [-] epsilon for numerical derivatives (finite difference perturbation)
        
        !> \var zero Threshold for zero (may be too high for some applications)
        !> Values |x| < zero treated as zero to avoid division by zero and underflow
        !> Default: 1×10⁻²⁰ [M/L³] - may need reduction for trace species studies
        !> Physical meaning: minimum detectable concentration (practical zero)
        !> Note: May need adjustment for ultra-trace species (e.g., radionuclides)
        real(kind=8) :: zero=1d-20 !< [M/L³] threshold below which values treated as zero
        
        !> \var control_factor Controls Delta_c1 in Newton algorithm
        !> Damping factor: x_new = x_old + control_factor * Δx
        !> Default: 0.1 [-] - strong damping for stability (10% of full Newton step)
        !> Physical meaning: fraction of Newton step to apply (prevents overshooting)
        !> Range: (0, 1] where 1 = full Newton, smaller values = more conservative
        !> Helps convergence for highly nonlinear systems (e.g., pH near buffering points)
        real(kind=8) :: control_factor=1d-1 !< [-] damping factor for Newton updates (0 < factor ≤ 1)
        
        !> \var niter_max Maximum number of iterations allowed
        !> Safety limit to prevent infinite loops in non-converging systems
        !> Default: 50 [#] - sufficient for most chemical systems
        !> Physical meaning: maximum attempts before declaring convergence failure
        !> Triggers error or time step reduction if exceeded
        integer(kind=4) :: niter_max=50 !< [#] maximum iterations before convergence failure
        
        !> \var k_div_max Maximum number of time step divisions allowed
        !> Limits time step refinement when convergence fails or stability violated
        !> Default: 5 [#] - allows Δt to be reduced to Δt/2^5 = Δt/32 for stiff problems
        !> Physical meaning: how many times can halve time step for stiff problems
        !> Prevents excessively small time steps (would make simulation impractically slow)
        integer(kind=4) :: k_div_max=5 !< [#] maximum time step subdivisions for stability
        
        !> \var est_prm Parameter to estimate kinetic reaction amounts for downstream waters
        !> Used in reactive transport to extrapolate kinetic reaction progress downstream
        !> Default: 0.0 [-] - no extrapolation (conservative)
        !> Physical meaning: weighting factor for estimating reaction amounts in advected water
        !> Helps initialize Newton solver for downstream cells using upstream solution
        real(kind=8) :: est_prm=0.0 !< [-] estimation parameter for downstream kinetic reactions (0 = no extrapolation)
    contains
        !> \brief Set estimation parameter for downstream kinetic reactions
        procedure :: read_CV_params
    end type CV_params_s

    contains
        !> \brief Reads convergence parameters from input file unit
        subroutine read_CV_params(this, path, root)
            class(CV_params_s) :: this
            character(len=*), intent(in) :: path !< Path to input files
            character(len=*), intent(in) :: root !< Root name for input files
            !integer(kind=4), intent(in) :: unit !< Input file unit number
            integer(kind=4) :: unit !< Input file unit number
            integer(kind=4) :: ierr !< I/O error flag
            character(len=256) :: line !< Full path to input file

            unit = 111 !< Arbitrary unit number for reading CV parameters
            open(unit,file=trim(path)//trim(root)//'_CV_params.dat',status='old',action='read',iostat=ierr)
            if (ierr /= 0) then
                write(*,*) 'Error: Could not open CV_params file: '//trim(path)//trim(root)//'_CV_params.dat'
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
                    ! read(unit,*) this%eps
                    ! read(unit,*) this%zero
                    ! read(unit,*) this%control_factor
                    ! read(unit,*) this%niter_max
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