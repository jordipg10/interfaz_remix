!> Computes variable activity species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
!! We assume all primary species are aqueous
!! The Jacobians are computed analytically
subroutine water_mixing_iter_EfI_kin_anal_ideal_lump(this,c1_old,c_tilde,conc_nc,porosity,Delta_t)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, initialise_iterative_method, &
        Newton_EfI_rk_kin_aq_anal_ideal
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:) !> aqueous variable activity concentrations at previous time step
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations (already allocated)
    real(kind=8), intent(in), optional :: porosity !> at this target
    real(kind=8), intent(in), optional :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: c1(:) !> aqueous var act concentrations
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    real(kind=8) :: Delta_t_bis !> reduced time step
    logical :: CV_flag !> convergence flag
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> We get aqueous var act species concentrations
        c1=this%get_conc_nc()
    !> Newton initialisation parameter
        mu=0d0 
!> Process
    !> Initial guess variabel activity concentrations
        call initialise_iterative_method(c1_old,c1,mu,conc_nc)
        call this%set_conc_var_act_species(conc_nc) !> chapuza
        call this%set_act_aq_species() !> chapuza
        k=0
        k_div=0
        Delta_t_bis=Delta_t
        !> Update old kinetic reaction rates
        call this%update_rk_old() !> we update old values of kinetic reaction rates
        call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
        do
            do !> loop until convergence is reached
            !> We apply Newton method to compute aqueous concentrations
                call this%Newton_EfI_rk_kin_aq_anal_ideal_lump(c_tilde,porosity,Delta_t_bis,conc_nc,niter,CV_flag)
            !> We check convergence
                if (CV_flag.eqv..false.) then !> no CV
                    if (mu<1d0) then
                        mu=mu+0.25 !> we increase Newton initialisation parameter
                        !call initialise_iterative_method(c1_old,c1,mu,conc_nc)
                        !call this%set_conc_var_act_species(conc_nc) !> chapuza
                        !call this%set_act_aq_species() !> chapuza
                    else
                        mu=0d0
                        k_div=k_div+1
                        Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
                    end if
                    call initialise_iterative_method(c1_old,c1,mu,conc_nc) !> chapuza
                    call this%set_conc_var_act_species(conc_nc) !> chapuza
                    call this%set_act_aq_species() !> chapuza
                else 
                    k=k+1
                    exit
                end if
            end do
            if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                exit
            end if
        end do
!< Post-process
    !call this%set_conc_var_act_species(conc_nc)
    deallocate(c1)
end subroutine