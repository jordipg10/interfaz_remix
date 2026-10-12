!> \file equilibrate_water_ideal.f90
!> \brief Speciate (equilibrate) a water with its reactive zone under the ideal-solution assumption.
!> \details Copy of mixing_iter_comp_ideal WITHOUT the "unchanged primaries" early-exit: the Newton
!> speciation ALWAYS runs. Used by the per-reactive-zone mixing (mix_within_zone), where the water
!> handed in is a displaced (transported) water whose current state is NOT yet at equilibrium with
!> the (new) reactive zone, so the early-exit of mixing_iter_comp_ideal would wrongly skip speciation.
subroutine equilibrate_water_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object (speciated in place)
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations from previous time step (Newton initial guess)
    real(kind=8), intent(in) :: c_hat(:) !> mixed variable activity species concentrations
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for kinetic reactions (unused; kept for signature parity)
    real(kind=8), intent(in) :: Delta_t !> time step (unused; kept for signature parity)
    real(kind=8), intent(in) :: theta !> time weighting factor (unused; kept for signature parity)
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations (output, already allocated)
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations (already allocated)

    integer(kind=4) :: niter !> Newton iterations performed
    logical :: CV_flag !> Newton convergence flag
    logical :: tried_c_hat !> .true. once the c_hat primaries have been tried as initial guess
    real(kind=8) :: mu !> Newton initialisation parameter in [0,1]
    real(kind=8), allocatable :: c1(:) !> this water's primaries on entry (retry extrapolation anchor; Newton overwrites this's primaries on failure)
    real(kind=8), allocatable :: c1_ig(:) !> Newton initial guess for primary species

    !> Component concentrations from the mixed species: u = U*c_hat
    conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_hat)
    !> ALWAYS speciate (no early-exit): Newton with adaptive initialisation mu = 0 -> 0.25 -> ... -> 1,
    !> then a last attempt from the primary block of c_hat itself (a conservative mix of equilibrated
    !> waters, so its primaries are usually close to the speciated ones).
    !> mu=0 gives c1_ig = c1 exactly, so the extrapolation is only evaluated on retries.
    !> compute_c_nc_from_u_Newton_ideal takes its guess from c1_ig (and sets this's primaries itself),
    !> so no set_conc_prim_species is needed here.
    c1=this%get_c1()
    mu=0d0
    tried_c_hat=.false.
    c1_ig=c1
    do
        call this%compute_c_nc_from_u_Newton_ideal(c1_ig,conc_comp,conc_nc,niter,CV_flag)
        if (CV_flag) exit
        if (mu<1d0) then
            mu=mu+0.25d0 !> retry extrapolating away from c1_old: c1_ig = (1+mu)*c1 - mu*c1_old
            call initialise_iterative_method(c1_old,c1,mu,c1_ig)
        else if (.not. tried_c_hat) then
            tried_c_hat=.true. !> last retry: the mixed primaries themselves
            c1_ig=c_hat(1:size(c1))
        else
            error stop "Newton speciation not converging, you must reduce time step"
        end if
    end do
end subroutine
