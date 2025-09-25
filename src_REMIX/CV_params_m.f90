!> Convergence parameters module
!! This structure contains numerical parameters related to convergence of iterative methods
module CV_params_m
    implicit none
    save
    type, public :: CV_params_s !> convergence parameters structure
        real(kind=8) :: abs_tol=1d-9 !> arithmetic absolute tolerance (residual Newton method)
        real(kind=8) :: log_abs_tol=1d-9 !> logarithmic absolute tolerance (residual Newton method)
        real(kind=8) :: rel_tol=1d-16 !> relative tolerance
        real(kind=8) :: eps=1d-12 !> epsilon for incremental coefficients
        real(kind=8) :: zero=1d-9 !> zero (demasiado alto)
        real(kind=8) :: control_factor=1d-1 !> controls Delta_c1 in Newton algorithm
        integer(kind=4) :: niter_max=80 !> maximum number of iterations
        integer(kind=4) :: k_div_max=6 !> maximum number of time step divisions
        real(kind=8) :: mu=0d0 !> parameter to estimate the kinetic reaction amounts of downstream waters
    end type
end module CV_params_m