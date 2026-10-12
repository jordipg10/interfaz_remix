!> Computes concentration of secondary species with variable activity from concentration of primary species using Picard method and mass action law
!! This subroutine is meant to be used when primary species are aqueous & solid
!!
!! Mass action law for the secondary variable activity species (log10):
!!   log c2v = Se_nc_1_star*(log gamma1 + log c1) + logK_star - log gamma2v
!! Only the activity coefficients of the AQUEOUS species change along the Picard iteration (they
!! depend on the ionic strength, i.e. on c2v itself); c1, logK_star and the activity coefficients of
!! the non-aqueous primary (solids, gases) and secondary (equilibrium var-act minerals, surface
!! complexes, gases) species are fixed, so their whole contribution is assembled ONCE before the loop
!! (log_c2v_cst). The iteration itself works on local arrays only (no object updates, allocations or
!! array temporaries per iteration):
!!   (1) ionic strength I = I_fixed + 1/2 sum_i z_i^2 c2v_i, where I_fixed (primary and constant
!!       activity aqueous species, from the object's concentrations) and the squared valences z_i^2 of
!!       the aqueous secondary variable activity species are precomputed;
!!   (2) log10 activity coefficients of the dissolved solids from I with the per-species model
!!       parameters gathered once: log gamma_i = -alpha_i A z_i^2 sqrt(I)/(1 + beta_i sqrt(I)) + gamma_i I
!!       (the formula of aq_phase_c%compute_log_act_coeffs_aq_phase; the water's log gamma is fixed);
!!   (3) the aqueous gamma corrections of the mass action law, accumulated column by column over the
!!       aqueous primary columns of Se_nc_1_star, and the exponentiation.
!! The object is updated once after convergence: ionic strength and log activity coefficients at the
!! converged c2v, the secondary concentrations, activities, the water activity coefficient and the
!! salinity (none of them is needed by the iteration itself).
!! If the fixed point has not converged after CV_params%niter_max iterations (e.g. the ionic strength
!! of an extreme trial point of a Newton line search makes the activity coefficients diverge) the
!! routine returns with CV_flag = .false. and leaves the object untouched; c2v then holds the last
!! (non-converged, possibly NaN) iterate. The callers decide how to react (the Newton methods reject
!! such a point as infeasible).
subroutine compute_c2v_from_c1_Picard(this,c1,c2v_ig,c2v,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1(:) !> primary concentrations (dim=n_p) (must be already set in the object)
    real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations (must be already set in the object)
    real(kind=8), intent(out) :: c2v(:) !> secondary variable activity concentrations (must be already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise (then the object is left untouched)
!> Variables
    integer(kind=4) :: n_p,n_p_aq,n_v2_aq,n_e,n_v_aq,n_aq,n_diss,niter_max,i,j,k,m
    real(kind=8) :: rel_tol,rel_change,c_zero,A_DH,I_fixed,ionic_str,lg_j
    real(kind=8), allocatable :: log_gamma1(:),log_gamma2v(:),log_c2v_cst(:),log_c2v_new(:),c2v_old(:)
    real(kind=8), allocatable :: lg(:) !> log10 activity coefficients of the aqueous species at the current ionic strength, concentration order (n_aq)
    real(kind=8), allocatable :: coef_A(:),coef_beta(:),coef_gamma(:) !> per dissolved solid: -alpha*A*z^2, beta, gamma (n_diss)
    integer(kind=4), allocatable :: ind_diss(:) !> concentration index of each dissolved solid (n_diss)
    integer(kind=4), allocatable :: ind_prim_aq(:) !> concentration index of each aqueous primary species (n_p_aq)
    integer(kind=4), allocatable :: ind_sec_aq(:) !> concentration index of each aqueous secondary variable activity species, c2v order (n_v2_aq)
    real(kind=8), allocatable :: z2_sec(:) !> squared valence of each aqueous secondary variable activity species (n_v2_aq)
    logical, allocatable :: is_sec(:) !> concentration index belongs to an aqueous secondary variable activity species (n_aq)
!> Pre-process
    associate(sa => this%solid_chemistry%reactive_zone%speciation_alg, Se => this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star, &
              aq => this%aq_phase)
    n_p=sa%num_prim_species
    n_p_aq=sa%num_aq_prim_species
    n_e=sa%num_eq_reactions
    n_v2_aq=sa%num_aq_sec_var_act_species
    n_v_aq=sa%num_aq_var_act_species
    n_aq=aq%num_species
    n_diss=n_aq-aq%wat_flag
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max
    c_zero=this%solid_chemistry%reactive_zone%CV_params%zero !> concentration treated as zero: floor of the relative-change denominator
    !> Constant part of the mass action law (fixed along the iteration):
    !>   Se_nc_1_star*log10(c1) + logK_star
    !>   + Se_nc_1_star(:,non-aqueous primaries)*log gamma1(non-aqueous primaries)
    !>   - log gamma2v(non-aqueous secondaries)
    log_gamma1=this%get_log_gamma1()
    log_gamma2v=this%get_log_gamma2v()
    allocate(log_c2v_cst(n_e),log_c2v_new(n_e),c2v_old(n_e))
    log_c2v_cst=matmul(Se,log10(c1))+this%logK_star
    if (n_p>n_p_aq) log_c2v_cst=log_c2v_cst+matmul(Se(:,n_p_aq+1:n_p),log_gamma1(n_p_aq+1:n_p))
    if (n_e>n_v2_aq) log_c2v_cst(n_v2_aq+1:n_e)=log_c2v_cst(n_v2_aq+1:n_e)-log_gamma2v(n_v2_aq+1:n_e)
    !> Concentration indices of the aqueous primary and secondary variable activity species (the latter
    !> in c2v order, as set_conc_sec_var_act_species scatters them) and the squared valences of the
    !> secondaries, the only ionic strength contributions that change along the iteration
    allocate(ind_prim_aq(n_p_aq),ind_sec_aq(n_v2_aq),z2_sec(n_v2_aq),is_sec(n_aq),lg(n_aq))
    ind_prim_aq=this%ind_prim_species(1:n_p_aq)
    ind_sec_aq=this%ind_var_act_species(n_p_aq+1:n_v_aq)
    is_sec=.false.
    do i=1,n_v2_aq
        m=ind_sec_aq(i)
        is_sec(m)=.true.
        z2_sec(i)=dble(aq%aq_species(this%indices_aq_phase(m))%valence**2)
    end do
    !> Fixed part of the ionic strength: primary and constant activity aqueous species (each
    !> concentration paired with the valence of ITS species through the water -> aqueous phase index map)
    I_fixed=0d0
    do m=1,n_aq
        if (is_sec(m)) cycle
        I_fixed=I_fixed+this%concentrations(m)*dble(aq%aq_species(this%indices_aq_phase(m))%valence**2)
    end do
    I_fixed=0.5d0*I_fixed
    !> Activity coefficient model parameters of the dissolved solids, gathered once (water excluded:
    !> its log gamma is kept as stored in the object, like compute_log_act_coeffs_aq_chem does).
    !> beta and gamma are evaluated with THIS water's Debye-Hückel A and B (temperature dependent).
    allocate(ind_diss(n_diss),coef_A(n_diss),coef_beta(n_diss),coef_gamma(n_diss))
    A_DH=this%params_aq_sol%A
    do i=1,n_diss
        k=aq%ind_diss_solids(i)
        ind_diss(i)=this%indices_aq_species(k)
        coef_A(i)=-(aq%aq_species(k)%params_act_coeff%alpha*(A_DH*dble(aq%aq_species(k)%valence**2)))
        coef_beta(i)=aq%aq_species(k)%params_act_coeff%get_beta(this%params_aq_sol)
        coef_gamma(i)=aq%aq_species(k)%params_act_coeff%get_gamma(this%params_aq_sol)
    end do
    lg=this%log_act_coeffs !> non-dissolved-solid entries (water) stay as stored
!> Process: Picard iteration on the aqueous activity coefficients
    c2v_old=c2v_ig
    ionic_str=I_fixed+0.5d0*dot_product(z2_sec,c2v_old) !> ionic strength of the initial guess
    call log_act_coeffs_diss_solids(ionic_str)
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    do
        niter=niter+1 !> we update number of iterations
    !> Mass action law: constant part + aqueous activity coefficient corrections (aqueous primary
    !> columns of Se_nc_1_star, then the aqueous secondaries' own gamma)
        log_c2v_new=log_c2v_cst
        do j=1,n_p_aq
            lg_j=lg(ind_prim_aq(j))
            if (lg_j/=0d0) log_c2v_new=log_c2v_new+Se(:,j)*lg_j
        end do
        do i=1,n_v2_aq
            log_c2v_new(i)=log_c2v_new(i)-lg(ind_sec_aq(i))
        end do
        c2v=1d1**log_c2v_new
    !> We check convergence (relative change of the secondary concentrations). The denominator is
    !> floored at the zero-concentration threshold so that a vanishing previous iterate (e.g. a zero
    !> initial guess) gives a finite change: 0 -> 0 counts as converged, 0 -> c>0 does not, instead of
    !> the Inf / NaN of a plain division. A NaN iterate (overflow) is never accepted as converged.
        rel_change=0d0
        do i=1,n_e
            if (c2v(i)/=c2v(i)) then
                rel_change=huge(1d0)
                exit
            end if
            rel_change=max(rel_change,abs(c2v(i)-c2v_old(i))/max(abs(c2v_old(i)),c_zero))
        end do
    !> Ionic strength and aqueous activity coefficients at the new iterate: next iteration's input,
    !> or the converged state stored in the object below
        ionic_str=I_fixed+0.5d0*dot_product(z2_sec,c2v)
        call log_act_coeffs_diss_solids(ionic_str)
        if (rel_change<rel_tol) then
            CV_flag=.true.
            exit !> CV reached
        else if (niter==niter_max) then
            exit !> no CV: leave with CV_flag = .false. (the caller decides; the object is not updated)
        end if
        c2v_old=c2v !> we update concentrations for next iteration
    end do
    end associate
!> Post-process: state consistent with the converged c2v (skipped when the fixed point failed, so that
!> the object keeps the state of the last converged speciation)
    if (CV_flag) then
        this%ionic_strength=ionic_str
        do i=1,n_diss
            this%log_act_coeffs(ind_diss(i))=lg(ind_diss(i)) !> log_10 activity coefficients of the dissolved solids at the converged ionic strength
        end do
        call this%set_conc_sec_var_act_species(c2v) !> secondary variable activity concentrations in the aqueous (and solid, gas) chemistry objects
        call this%compute_activities() !> we compute activities
        call this%compute_pH()
        call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
        call this%compute_salinity()
    end if
    deallocate(log_gamma1,log_gamma2v,log_c2v_cst,log_c2v_new,c2v_old,lg,coef_A,coef_beta,coef_gamma,ind_diss,ind_prim_aq, &
        ind_sec_aq,z2_sec,is_sec)
contains
    !> log10 activity coefficients of the dissolved solids at ionic strength I (stored in lg by
    !> concentration index): -alpha*A*z^2*sqrt(I)/(1+beta*sqrt(I)) + gamma*I, as compute_log_act_coeffs_aq_phase
    subroutine log_act_coeffs_diss_solids(I)
        real(kind=8), intent(in) :: I
        real(kind=8) :: sqrt_I
        integer(kind=4) :: ii
        sqrt_I=sqrt(I)
        do ii=1,n_diss
            lg(ind_diss(ii))=coef_A(ii)*sqrt_I/(1d0+coef_beta(ii)*sqrt_I)+coef_gamma(ii)*I
        end do
    end subroutine log_act_coeffs_diss_solids
end subroutine compute_c2v_from_c1_Picard