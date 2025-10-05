!> Computes component and variable activity concentrations after iteration of WMA method for equilibrium chemical system
!> We assume all primary species are aqueous & solid
subroutine mixing_iter_comp_exch_ideal(this,c1_old,c_tilde,conc_nc,porosity,Delta_t)
    use aqueous_chemistry_m,    only: aqueous_chemistry_c, inf_norm_vec_real, initialise_iterative_method
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    real(kind=8), intent(in) :: c_tilde(:)
    real(kind=8), intent(out) :: conc_nc(:)
    !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
    real(kind=8), intent(in), optional :: porosity
    real(kind=8), intent(in), optional :: Delta_t !> time step
!> Variables
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8), allocatable :: c1(:) !> concentration primary species
    real(kind=8), allocatable :: c1_ig(:) !> initial guess concentration primary species
    real(kind=8), allocatable :: conc_comp(:) !> concentration components
!> Pre-process
    c1=this%get_c1()
    allocate(c1_ig(THIS%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
!> Process  
    !> We compute component concentrations after mixing
        conc_comp=MATMUL(THIS%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_tilde)
    !> Loop until speciation converges
        do
        !> We initialise primary concentrations for Newton speciation
            call initialise_iterative_method(c1_old,c1,mu,c1_ig)
        !> We compute variable activity concentrations from component concentrations
            call this%compute_c_nc_from_u_Newton_ideal(c1_ig,conc_comp,conc_nc,niter,CV_flag)
        !> We check convergence
            if (CV_flag .eqv. .false.) then !> NO CV
                if (mu<1d0) then
                    mu=mu+0.25
                else
                    error stop "Newton speciation does not converge"
                end if
            else
                exit
            end if
        end do
        call this%solid_chemistry%compute_activities_solids()
!> Post-process
    !> We check concentrations
        !call this%check_conc_aq_var_act_species(conc_comp)
        !call this%check_act_aq_species()
end subroutine