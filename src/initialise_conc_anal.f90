!> \file initialise_conc_anal.f90
!> \brief Newton speciation with analytical Jacobian (general activity coefficients).
!> \details
!> Computes species concentrations from water-type constraints using
!> Newton's method with the analytical Jacobian, accounting for
!> non-ideal activity coefficients.  Used during CHEPROO water-type reading.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup chemistry
!> \see initialise_conc_anal_ideal, initialise_conc_incr_coeff, aqueous_chemistry_m

!> Computes species concentrations with data from water type definition using Newton method
!> We assume the initial guess of primary aqueous species is already set in the aqueous chemistry object
!> This subroutine is only used to read water types based in CHEPROO
subroutine initialise_conc_anal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use arrays_m, only: int_array_c
    use metodos_sist_lin_m, only: LU_lin_syst
    use vectors_m, only: inf_norm_vec_real, outer_prod_vec
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: icon(:) !> initial condition type
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> indices of constrains in equilibrium stoichiometric matrix
    real(kind=8), intent(in) :: ctot(:) !> data given
    integer(kind=4), intent(out) :: niter !> number of iterations Newton-Raphson
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    !character(len=*), intent(in), optional :: path_py !> path output Python
!> Variables
    real(kind=8), allocatable :: c2_old(:),c2_new(:),conc_comp(:),c1(:),log_c2k(:),log_c2(:),dc2_dc1(:,:),out_prod(:,:),&
        log_gamma_ideal(:)
    real(kind=8), allocatable :: res(:) !> residual in Newton-Raphson
    real(kind=8), allocatable :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: tol_res(:) !> tolerance residues Newton-Raphson
    real(kind=8), allocatable :: mat_lin_syst(:,:),Se_aq_comp(:,:),K(:),u_aq(:),z2(:),d_log_gamma_d_I(:),&
        log_Jacobian_act_coeffs(:,:)
    real(kind=8), allocatable :: z2_sp(:) !> squared valences in speciation-algebra (species) order (0 for non-aqueous species)
    real(kind=8), allocatable :: c1_save(:),scaled_Delta(:),res_trial(:),Jac_trial(:,:) !> line search: saved primaries, scaled step, trial residual, discarded trial Jacobian
    real(kind=8) :: res_norm,alpha,scale_factor !> line search: current residual norm, step length, step bound factor
    integer(kind=4) :: i_LS !> line search trial counter
    integer(kind=4), parameter :: max_LS_tries=10 !> maximum number of step halvings
    real(kind=8), parameter :: max_dp=2d0 !> maximum change per Newton step, in decades (ln(10)*max_dp*c bounds |Delta_c|)
    logical :: step_accepted
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    integer(kind=4) :: niter_ideal !> Newton iterations of the ideal start (fallback)
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:)
    logical :: CV_flag_Picard,zero_flag
    logical :: CV_flag_ideal !> convergence of the ideal start (fallback)
    logical :: ideal_start_done !> the ideal start has been tried (only once)
    type(int_array_c) :: indices_icon
    
!> Pre-process
    CV_flag=.false.
    ideal_start_done=.false.
    if (size(icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) error stop
    
    if (sum(n_icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) error stop
    
    allocate(c2_old(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),&
    c2_new(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    ind_cstr=0 !> index of phase constrain
    call indices_icon%allocate_array(4) !> number of icon options
    do i=1,indices_icon%num_cols
        call indices_icon%cols(i)%allocate_vector(n_icon(i))
    end do
    allocate(counters(4)) !> dim=n�icon options
    counters=0 !> initial counter for each icon option
    allocate(tol_res(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        tol_res(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol !> absolute tolerance (by default)
        !> Concentrations
        if (icon(i)==1) then
            counters(1)=counters(1)+1
            indices_icon%cols(1)%col_1(counters(1))=i
        !> Aqueous components
        else if (icon(i)==2) then
            counters(2)=counters(2)+1
            indices_icon%cols(2)%col_1(counters(2))=i
        !> Activities
        else if (icon(i)==3) then
            counters(3)=counters(3)+1
            indices_icon%cols(3)%col_1(counters(3))=i
        !> Phases
        else if (icon(i)==4) then
            counters(4)=counters(4)+1
            indices_icon%cols(4)%col_1(counters(4))=i
            tol_res(i)=1d1**this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !> phase constraint: the residual is logarithmic (sum Se log10 a + log10 ctot - logK), so its tolerance is 10**log_rel_tol (relative error ~ ln(10)*tol on the constrained activity), NOT log_rel_tol itself (a negative number whose magnitude made the test vacuous)
        else
            error stop "icon option not implemented"
        end if
    end do
    !> tol_res was built by INPUT position, i.e. by the concentration index k of the constrained primary
    !> (icons(k), ctots(k)); res and Jac_res (compute_res_Jac_res_anal) come out in the unknowns'
    !> stoichiometric order, row j = primary stored at ind_prim_species(j). Permute the tolerances to the
    !> rows so that each residual is tested against ITS tolerance: a water type that lists its species in
    !> another order than the chemical system would otherwise test a logarithmic phase residual against
    !> abs_tol (which it can never meet) and stop with "Newton method not accurate enough".
    tol_res=tol_res(this%ind_prim_species(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    !print *, indices_icon%cols(4)%col_1
!> Newton-Raphson
    allocate(res(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),&
    Jac_res(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    allocate(Delta_c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
    c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_species),&
        z2_sp(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
    !> rows and columns indexed as compute_res_Jac_res_anal reads them (concentration index / shifted secondaries), see below
    allocate(log_Jacobian_act_coeffs(this%aq_phase%num_species,max(this%aq_phase%num_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)))
    allocate(dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> chapuza    
    allocate(c1_save(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
        scaled_Delta(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
        res_trial(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),&
        Jac_trial(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    
    c2_old=1d-16 !> chapuza
    
    niter=0
    !this%log_act_coeffs=0d0 !> chapuza
    d_log_gamma_d_I=0d0 !> chapuza
    log_Jacobian_act_coeffs=0d0 !> mchapuza
    log_gamma_ideal=this%compute_log_act_coeffs_ideal()
    !> Squared charges species
    !allocate(z2(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
    !z2=0d0
    !z2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)=this%aq_phase%z2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
    !z2(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)=this%aq_phase%z2(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)
    do
    !> We update number of iterations
        niter=niter+1
        if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            print *, '[DEBUG initialise_conc_anal] Newton failed after', niter-1, 'iters'
            print *, '  ||res||_inf =', inf_norm_vec_real(res)
            print *, '  c1 =', c1
            error stop "Too many Newton iterations in initialisation"
        end if
        !call this%compute_ionic_strength()
        !call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_strength,this%params_aq_sol,this%log_act_coeffs)
        c1=this%get_c1()
        call this%compute_c2_from_c1_ideal(c1,log_gamma_ideal,c2_old)
        call this%compute_c2_from_c1_Picard(c1,c2_old,c2_new,niter_Picard,CV_flag_Picard)
        if (.not. CV_flag_Picard) then
        !> The activity coefficient fixed point does not exist at this iterate: with a poor input guess
        !> (e.g. h+ = 1e-9 against 0.25 M acetic acid: acetate = K*c(acid)/c(h+) ~ 1e6 M) the ionic
        !> strength makes the Davies/Debye-Hückel terms overflow. The iterate is only unchecked at the
        !> first iteration (the line search below rejects infeasible trial points), so recover ONCE by
        !> taking the IDEAL solution (gamma = 1, well behaved from any guess) as the starting point of
        !> the non-ideal Newton; a second failure is a genuine one.
            if (niter==1 .and. .not. ideal_start_done) then
                ideal_start_done=.true.
                write(0,'(3A)') ' NOTE: initialise_conc_anal (', trim(this%name), &
                    '): the activity coefficient fixed point diverges at the input guess; starting from the ideal (gamma = 1) solution'
                call this%initialise_conc_anal_ideal(icon,n_icon,indices_constrains,ctot,niter_ideal,CV_flag_ideal)
                if (.not. CV_flag_ideal) error stop "initialise_conc_anal: the ideal start did not converge"
                niter=0 !> the non-ideal iteration count restarts from the ideal solution
                cycle
            end if
            print *, '[initialise_conc_anal] Picard speciation did not converge at a Newton iterate (', trim(this%name), ')'
            print *, '  c1 =', c1
            error stop "Picard speciation failed in initialisation"
        end if
        !call this%compute_ionic_strength() !> we compute ionic activity
        !call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_strength,this%params_aq_sol,this%log_act_coeffs) !> we compute log activity coefficients aqueous species
        !call this%compute_activities()
        !call this%compute_log_act_coeff_wat()
        
        !print *, 10**this%log_act_coeffs
        !print *, c2_new
        !print *, this%activities
    !> Chapuza
        !c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)=this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)
    !> We compute residue and Jacobian analytically
        !> d(log gamma)/dI of every species in the SPECIATION-ALGEBRA order (the n_p primaries, then the
        !> n_e secondaries = the rows of Se_1_star, i.e. the order of c1 and c2): compute_dc2_dc1 reads
        !> out_prod in that order. compute_d_log_gamma_d_I_aq_chem is not used: it returns only the
        !> aqueous variable activity species and reads the aqueous phase BY POSITION, which is only right
        !> when this water's species are in aqueous-phase order (not after a primary swap); here the
        !> derivative is evaluated per species through the water's index maps. Non-aqueous species
        !> (solids, gases) and water have no ionic-strength dependence: 0.
        !> The outer product with z^2 (also species-ordered) gives d(log gamma_i)/dI * z_j^2.
        d_log_gamma_d_I=0d0
        z2_sp=0d0
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            call d_log_gamma_d_I_species(this%ind_prim_species(j),d_log_gamma_d_I(j),z2_sp(j))
        end do
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            call d_log_gamma_d_I_species(this%ind_sec_species(j),&
                d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j),&
                z2_sp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j))
        end do
        out_prod=outer_prod_vec(d_log_gamma_d_I,z2_sp)
        !> We compute Jacobian secondary-primary concentrations
        call this%compute_dc2_dc1(out_prod,c1,c2_new,dc2_dc1)
        !> Log-Jacobian activity coefficients-concentrations, d(log gamma_i)/d(log c_j) = 0.5 ln(10) c_j
        !> d(log gamma_i)/dI z_j^2, in the layout compute_res_Jac_res_anal reads it: rows and primary
        !> columns by CONCENTRATION index, aqueous secondary columns at num_prim_species+1:num_aq_species
        !> (paired with the rows of dc2aq_dc1). Only aqueous species have nonzero entries.
        log_Jacobian_act_coeffs=0d0
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_Jacobian_act_coeffs(this%ind_prim_species(i),this%ind_prim_species(j))=5d-1*log(1d1)*&
                    this%concentrations(this%ind_prim_species(j))*out_prod(i,j)
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                log_Jacobian_act_coeffs(this%ind_sec_species(i),this%ind_prim_species(j))=5d-1*log(1d1)*&
                    this%concentrations(this%ind_prim_species(j))*&
                    out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i,j)
            end do
        end do
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                log_Jacobian_act_coeffs(this%ind_prim_species(i),&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j)=5d-1*log(1d1)*&
                    this%concentrations(this%ind_sec_species(j))*&
                    out_prod(i,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j)
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                log_Jacobian_act_coeffs(this%ind_sec_species(i),&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j)=5d-1*log(1d1)*&
                    this%concentrations(this%ind_sec_species(j))*&
                    out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i,&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+j)
            end do
        end do
    !> We check dc2_dc1
        !> NOTE: check_dc2_dc1 is a development-time validation of the analytical
        !> Jacobian. It is disabled here because (a) its threshold uses log_rel_tol
        !> directly (a negative log-scale value) against a non-negative matrix norm,
        !> so it can never pass, and (b) it is an absolute comparison that is not
        !> meaningful for non-ideal systems where the Jacobian entries can be large
        !> away from the solution. The Newton residual test below is the actual
        !> convergence criterion.
        !call this%check_dc2_dc1(c1,c2_new,dc2_dc1,log_Jacobian_act_coeffs)
        call this%compute_res_Jac_res_anal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species,:),log_Jacobian_act_coeffs,res,Jac_res)
    !> We check convergence: every constraint residual within its own tolerance (concentration residuals
    !> against abs_tol, logarithmic phase residuals against 10**log_rel_tol)
        if (all(abs(res)<=tol_res)) then !> CV reached
            CV_flag=.true.
            exit
        end if
    !> We solve linear system aqueous primary concentartions
        call LU_lin_syst(Jac_res,-res,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1)
        !> c1^(i+1)=c1^i+Delta_c1^i
        if (inf_norm_vec_real(Delta_c1/c1)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            print *, '[initialise_conc_anal] Newton step below rel_tol while the residual is still above its tolerance (', &
                trim(this%name), ')'
            print *, '  ||Delta_c1/c1||_inf =', inf_norm_vec_real(Delta_c1/c1)
            print *, '  res     =', res, ' (unknowns in stoichiometric order)'
            print *, '  tol_res =', tol_res
            print *, '  c1      =', c1
            error stop "Newton method not accurate enough in initialisation"
        end if
    !> Globalisation (as initialise_conc_anal_ideal): the step is bounded to max_dp decades per species
    !> (uniform scaling, direction preserved) and its length alpha is halved until the residual norm
    !> decreases (backtracking on ||res||_inf; the trial residual needs the Picard speciation at the
    !> trial point, evaluated by residual_at). A raw Newton step from a poor water-type guess (e.g.
    !> h+ = 1e-7 against a calcite constraint) otherwise oscillates between c and c/10 for ever.
        res_norm=inf_norm_vec_real(res)
        c1_save=c1
        scale_factor=1d0
        do i=1,size(Delta_c1)
            if (abs(Delta_c1(i))>0d0 .and. c1_save(i)>0d0) scale_factor=min(scale_factor,max_dp*log(1d1)*c1_save(i)/abs(Delta_c1(i)))
        end do
        if (scale_factor<1d0) Delta_c1=Delta_c1*scale_factor
        alpha=1d0
        step_accepted=.false.
        do i_LS=1,max_LS_tries
            scaled_Delta=alpha*Delta_c1
            c1=c1_save
            call this%update_conc_aq_prim_species(c1,scaled_Delta,zero_flag) !> c1 <- c1_save + alpha*Delta (positivity guard inside)
            call residual_at(c1,res_trial)
            if (inf_norm_vec_real(res_trial)<res_norm) then
                step_accepted=.true.
                exit
            end if
            alpha=alpha*0.5d0
        end do
        if (.not. step_accepted) then
            print *, '[initialise_conc_anal] no step reduces the residual; ||res||_inf =', res_norm
            print *, '  c1 =', c1_save
            print *, '  res =', res
            error stop "Newton line search failed in initialisation"
        end if
        c2_old=c2_new
    end do
    ! call this%compute_pH()
    ! call this%compute_salinity()
    ! call this%compute_alkalinity()
    !if (associated(this%solid_chemistry)) then
    !    nullify(this%solid_chemistry)
    !end if
    !if (associated(this%gas_chemistry)) then
    !    nullify(this%gas_chemistry)
    !end if
contains
    !> Residual of the water-type constraints at the primary concentrations c1_try (already set in the
    !> object by update_conc_aq_prim_species): ideal guess -> Picard speciation -> compute_res_Jac_res_anal.
    !> Only the residual is used (the Jacobian inputs are the current iterate's; the Jacobian output is
    !> discarded). Leaves the object speciated at c1_try, which is the next iterate when the step is
    !> accepted; a rejected trial is overwritten by the next trial / the accepted one.
    subroutine residual_at(c1_try,res_try)
        real(kind=8), intent(in) :: c1_try(:)
        real(kind=8), intent(out) :: res_try(:)
        call this%compute_c2_from_c1_ideal(c1_try,log_gamma_ideal,c2_old)
        call this%compute_c2_from_c1_Picard(c1_try,c2_old,c2_new,niter_Picard,CV_flag_Picard)
        if (.not. CV_flag_Picard) then
            res_try=huge(1d0) !> infeasible trial point (activity coefficient fixed point diverged): never accepted by the line search
            return
        end if
        call this%compute_res_Jac_res_anal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species,:),log_Jacobian_act_coeffs,res_try,Jac_trial)
    end subroutine residual_at

    !> d(log10 gamma)/dI and z^2 of the aqueous species stored at concentration index k, at the current
    !> ionic strength: derivative of the model of compute_log_act_coeffs_aq_phase,
    !> log gamma = -alpha*A*z^2*sqrt(I)/(1+beta*sqrt(I)) + gamma*I, i.e.
    !> d(log gamma)/dI = -alpha*A*z^2/(2*sqrt(I)*(1+beta*sqrt(I))^2) + gamma (as compute_d_log_gamma_d_I_aq_chem).
    !> Water (constant activity in the model) gets 0. The species is reached through this water's
    !> concentration -> aqueous phase index map, so a reordered (swapped) water is handled correctly.
    subroutine d_log_gamma_d_I_species(k,dlg,z2k)
        integer(kind=4), intent(in) :: k !> concentration index of the species
        real(kind=8), intent(out) :: dlg !> d(log10 gamma)/dI
        real(kind=8), intent(out) :: z2k !> squared valence
        integer(kind=4) :: ka !> index in the aqueous phase
        real(kind=8) :: sqrt_I
        ka=this%indices_aq_phase(k)
        z2k=dble(this%aq_phase%aq_species(ka)%valence**2)
        if (ka==this%aq_phase%ind_wat) then
            dlg=0d0
            return
        end if
        sqrt_I=sqrt(this%ionic_strength)
        associate(pac => this%aq_phase%aq_species(ka)%params_act_coeff)
        dlg=-pac%alpha*(this%params_aq_sol%A*z2k)/(2d0*sqrt_I*(1d0+pac%get_beta(this%params_aq_sol)*sqrt_I)**2)+ &
            pac%get_gamma(this%params_aq_sol)
        end associate
    end subroutine d_log_gamma_d_I_species
end subroutine