!> Computes aqueous species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
!! We assume all primary species are aqueous
!! The Jacobians are computed analytically
!! We apply lumping to the mixing ratios of kinetic reaction rates
subroutine water_mixing_iter_EfI_eq_kin_anal_ideal_lump(this,c1_old,c_tilde,Delta_t,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, initialise_c1_aq_iterative_method, &
        Newton_EfI_rk_eq_kin_aq_anal_ideal_lump
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:) !> aqueous primary concentrations at previous time step
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations (already allocated)
    !real(kind=8), intent(in), optional :: porosity !> at this target
!> Variables
    real(kind=8), allocatable :: c1(:) !> aqueous primary concentrations
    real(kind=8), allocatable :: u_tilde(:) !> component concentrations after mixing
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    real(kind=8) :: Delta_t_bis !> reduced time step
    logical :: CV_flag !> convergence flag
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> We get aqueous primary species concentrations
        c1=this%get_c1()
    !> Newton initialisation parameter
        mu=0d0 
!> Process
    !> Initial guess aquoeus primary concentrations
        call this%initialise_c1_aq_iterative_method(mu)
        k=0 !> counter of time steps
        k_div=0 !> counter of time step divisions
        Delta_t_bis=Delta_t
        !call this%compute_molalities() !> we change units to compute activities
        !> Update old kinetic reaction rates (chapuza)
        call this%update_rk_old() !> we update old values of kinetic reaction rates
        call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
     !> We compute component concentrations after mixing
        u_tilde=this%compute_u_tilde(c_tilde)
        do
            do !> loop until convergence is reached
            !> We apply Newton method to compute aqueous concentrations
                call this%Newton_EfI_rk_eq_kin_aq_anal_ideal_lump(u_tilde,Delta_t_bis,conc_nc,niter,CV_flag)
            !> We check convergence
                if (CV_flag.eqv..false.) then !> no CV
                    if (mu<1d0) then
                        mu=mu+0.25 !> we increase Newton initialisation parameter
                    else
                        mu=0d0
                        k_div=k_div+1
                        if (k_div>this%solid_chemistry%reactive_zone%CV_params%k_div_max) error stop "Too many time step divisions"
                        Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
                    end if
                    call this%initialise_c1_aq_iterative_method(mu)
                else 
                    k=k+1
                    exit
                end if
            end do
            if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                exit
            end if
        end do
    !> We compute component concentrations
        !conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
!< Post-process
        deallocate(c1,u_tilde)
end subroutine