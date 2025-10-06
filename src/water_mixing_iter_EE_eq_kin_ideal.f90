!> Computes aqueous variable activity species concentrations and aqueous component concentrations after iteration of WMA-EE method for equilibrium-kinetic chemical system
!! We assume all primary species are aqueous
subroutine water_mixing_iter_EE_eq_kin_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method, reaction_iteration_EE_eq_kin, &
    compute_u_tilde, compute_c_nc_from_u_aq_Newton_ideal, check_conc_var_act_species, check_act_aq_species
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations in previous time step
    !real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
    real(kind=8), intent(in) :: mix_ratio_Rk !> NO SE USA
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations after time iteration
    !real(kind=8), intent(in), optional :: porosity !> porosity in this target
    !real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: u_tilde(:)
    real(kind=8), allocatable :: conc_comp_react(:)
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: c1(:)
    real(kind=8), allocatable :: conc_nc_old(:)
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: k_div !> counter of time step divisions
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8) :: Delta_t_aux !> time step (auxiliary)
!> Pre-process
    allocate(conc_comp_react(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(conc_comp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(conc_nc_old(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    c1=this%get_c1() !> primary species concentrations
    conc_nc_old(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=c1 !> chapuza
    Delta_t_aux=Delta_t !> we set default value for time step
    k=0 !> we initialise counter of time steps
    k_div=0 !> we initialise counter of time step divisions
!> Process
    !> Components concentrations
    u_tilde=this%compute_u_tilde(c_tilde) !> transport part: concentrations of variable activity species after mixing
    conc_comp=u_tilde
    ! !> Update old kinetic reaction rates
    ! call this%update_rk_old() !> we update old values of kinetic reaction rates
    ! call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> We initialise primary concentrations for Newton speciation
    call this%initialise_c1_aq_iterative_method(mu)
    do
        !k=0
        !k_div=0
        !Delta_t_bis=Delta_t
        call this%reaction_iteration_EE_eq_kin(Delta_t_aux,rk_tilde,conc_comp_react) !> chemical part
        conc_comp=conc_comp+conc_comp_react !> we sum both parts
    !> Speciation for variable activity species concentrations
        !> Loop until speciation converges
        do
        !> We compute variable activity concentrations from component concentrations
            call this%compute_c_nc_from_u_aq_Newton_ideal(conc_comp,conc_nc,niter,CV_flag)
        !> We check convergence
            if (.not. CV_flag) then !> NO CV
                if (mu<1d0) then
                    mu=mu+0.25
                    call this%initialise_c1_aq_iterative_method(mu)
                else
                    mu=0d0
                    k_div=k_div+1
                    if (k_div>this%solid_chemistry%reactive_zone%CV_params%k_div_max) error stop "Too many time step divisions"
                    Delta_t_aux=Delta_t_aux/2d0 !> we reduce time step
                    conc_comp=conc_comp-conc_comp_react !> we reset component concentrations
                    call this%initialise_c1_aq_iterative_method(mu)
                    exit
                    !error stop "Newton speciation does not converge"
                end if
                
            else
                conc_nc_old=conc_nc
                k=k+1
                exit
            end if
        end do
        if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
            exit
        end if
    end do
        !> We check results
            !call this%check_conc_var_act_species(conc_nc,conc_comp)
            !call this%check_act_aq_species()
!> Post-process
    deallocate(conc_comp_react,conc_comp,c1,u_tilde)
end subroutine