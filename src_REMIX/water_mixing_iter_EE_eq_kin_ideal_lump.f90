!> Computes aqueous variable activity species concentrations and aqueous component concentrations after iteration of WMA-EE method for equilibrium-kinetic chemical system
!! We assume all primary species are aqueous
!! We apply lumping
subroutine water_mixing_iter_EE_eq_kin_ideal_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method, reaction_iteration_EE_eq_kin, &
    compute_u_tilde, compute_c_nc_from_u_aq_Newton_ideal, check_conc_var_act_species, check_act_aq_species
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations in previous time step
    !real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
    !real(kind=8), intent(in) :: mix_ratio_Rk !> NO SE USA
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations after time iteration
    !real(kind=8), intent(in), optional :: porosity !> porosity in this target
    
!> Variables
    real(kind=8), allocatable :: u_tilde(:)
    real(kind=8), allocatable :: conc_comp_react(:)
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: c1(:)
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
!> Pre-process
    allocate(conc_comp_react(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(conc_comp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    c1=this%get_c1() !> primary species concentrations
!> Components concentrations
    u_tilde=this%compute_u_tilde(c_tilde) !> transport part: concentrations of variable activity species after mixing
    !> Update old kinetic reaction rates
    !call this%update_rk_old() !> we update old values of kinetic reaction rates
    !call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
    
    call this%reaction_iteration_EE_eq_kin_lump(Delta_t,conc_comp_react) !> chemical part
    conc_comp=u_tilde+conc_comp_react !> we sum both parts
    
!> Speciation for variable activity species concentrations
    !> Loop until speciation converges
        do
        !> We initialise primary concentrations for Newton speciation
            call this%initialise_c1_aq_iterative_method(mu)
        !> We compute variable activity concentrations from component concentrations
            call this%compute_c_nc_from_u_aq_Newton_ideal(conc_comp,conc_nc,niter,CV_flag)
        !> We check convergence
            if (.not. CV_flag) then !> NO CV
                if (mu<1d0) then
                    mu=mu+0.25
                else
                    error stop "Newton speciation does not converge"
                end if
            else
                exit
            end if
        end do
    !> We check results
        !call this%check_conc_var_act_species(conc_nc,conc_comp)
        !call this%check_act_aq_species()
!> Post-process
    deallocate(conc_comp_react,conc_comp,c1,u_tilde)
end subroutine