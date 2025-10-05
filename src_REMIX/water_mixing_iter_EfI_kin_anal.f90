!> Computes aqueous species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
!! We assume all species are aqueous
!! The Jacobians are computed analytically
subroutine water_mixing_iter_EfI_kin_anal(this,c1_old,c2nc_ig,c_tilde,conc_nc,conc_comp,porosity,Delta_t)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:)
    real(kind=8), intent(out) :: conc_nc(:)
    real(kind=8), intent(out) :: conc_comp(:)
    real(kind=8), intent(in), optional :: porosity
    real(kind=8), intent(in), optional :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: conc(:) !> aqueous concentrations
    integer(kind=4) :: n_aq !> number of aqueous species
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    logical :: CV_flag !> convergence flag
!> Pre-process
        n_aq=this%aq_phase%num_species
    !> We assign aqueous species concentrations
        conc=this%concentrations
    !> Newton initialisation parameter
        mu=0d0
    !> Update old kinetic reaction rates
    call this%update_rk_old() !> we update old values of kinetic reaction rates
    call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> Process
    !> Loop until convergence is reached
        do
        !> We initialise aqueous concentrations
            call initialise_iterative_method([c1_old,c2nc_ig],conc,mu,this%concentrations)
        !> We apply Newton method to compute aqueous concentrations
            !call this%Newton_EfI_rk_kin_aq_anal(c_tilde,porosity,Delta_t,niter,CV_flag)
        !> We check convergence
            if (.not. CV_flag) then !> no CV
                if (mu<1d0) then
                    mu=mu+0.25 !> we increase Newton initialisation parameter
                else
                    error stop "Iterative method does not converge"
                end if
            else 
                exit
            end if
        end do
    conc_nc=this%concentrations
    conc_comp=this%concentrations
!> Post-process
    deallocate(conc)
end subroutine