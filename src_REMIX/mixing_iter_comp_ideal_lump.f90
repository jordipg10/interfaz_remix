!> Computes component and variable activity concentrations after iteration of WMA method for equilibrium chemical system
!> We assume all primary species are aqueous
subroutine mixing_iter_comp_ideal_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
    use aqueous_chemistry_m,   only: aqueous_chemistry_c, inf_norm_vec_real
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:) !> primary concentrations at previous time step
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
    real(kind=8), intent(out) :: conc_nc(:)
    
!> Variables
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8) :: Delta_t_bis !> reduced time step
    real(kind=8), allocatable :: c1(:) !> concentration primary species
    real(kind=8), allocatable :: conc_comp(:) !> concentration components
!> Pre-process
    c1=this%get_c1()
!> Process
    !> We compute component concentrations after mixing
        conc_comp=MATMUL(THIS%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_tilde)
        if (inf_norm_vec_real(c_tilde(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)-c1)<&
        this%solid_chemistry%reactive_zone%CV_params%zero) then
            conc_nc=this%get_conc_nc()
        else
        !> We initialise primary concentrations for Newton speciation
            call this%initialise_c1_aq_iterative_method(mu)
            !k=0
            !k_div=0
            !Delta_t_bis=Delta_t
        !> Loop until speciation converges
            !do
                do
                !> We compute variable activity concentrations from component concentrations
                    call this%compute_c_nc_from_u_aq_Newton_ideal(conc_comp,conc_nc,niter,CV_flag)
                !> We check convergence
                    if (CV_flag .eqv. .false.) then !> no CV
                        if (mu<1d0) then
                            mu=mu+0.25 !> we increase Newton initialisation parameter
                            call this%initialise_c1_aq_iterative_method(mu)
                        else
                            !mu=0d0
                            !k_div=k_div+1
                            !Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
                            error stop "Newton speciation not converging, you must reduce time step"
                        end if
                    else 
                        k=k+1
                        exit
                    end if
                end do
                !if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                !    exit
                !end if
            !end do
        end if
!> Post-process
    !> We check concentrations
        !call this%check_conc_aq_var_act_species(conc_comp)
        !call this%check_act_aq_species()
end subroutine