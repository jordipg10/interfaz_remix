!> Computes component and variable activity concentrations after iteration of WMA method for equilibrium chemical system
!> We assume primary species are aqueous & solid
subroutine mixing_iter_comp_exch(this,c1_old,c2nc_ig,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
    use aqueous_chemistry_m,    only: aqueous_chemistry_c, inf_norm_vec_real, initialise_iterative_method
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:)
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing, ordered in upstream & downstream respectively
    real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_nc(:)
    !real(kind=8), intent(in), optional :: porosity
    !real(kind=8), intent(in), optional :: Delta_t !> time step
!> Variables
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    integer(kind=4) :: n_p !> number of iterations in Newton speciation
    integer(kind=4) :: k !> number of iterations in Newton speciation
    integer(kind=4) :: k_div !> number of iterations in Newton speciation
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8) :: Delta_t_bis !> Newton initialistaion parameter
    real(kind=8), allocatable :: c1(:) !> concentration primary species
    real(kind=8), allocatable :: c1_ig(:) !> initial guess concentration primary species
    real(kind=8), allocatable :: conc_comp(:) !> concentration components
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    allocate(c1_ig(n_p))
    !> We assign primary species concentrations
        c1=this%get_c1()
    !> Newton initialisation parameter     
        mu=0d0 
!> Process
    !> We compute component concentrations after mixing
        conc_comp=MATMUL(THIS%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_tilde)
        if (inf_norm_vec_real(c_tilde(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)-c1(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))<&
        this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> chapuza
            conc_nc=this%get_conc_nc()
        else
            !> Loop until speciation converges
            do
            !> We initialise primary concentrations for Newton speciation
                call initialise_iterative_method(c1_old,c1,mu,c1_ig)
            !> We compute variable activity concentrations from component concentrations
                call this%compute_c_nc_from_u_Newton(c1_ig,c2nc_ig,conc_comp,conc_nc,niter,CV_flag)
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
        end if
!< Post-process
    deallocate(c1,c1_ig,conc_comp)
end subroutine