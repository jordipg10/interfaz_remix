!> Computes aqueous species concentrations from component concentrations using Newton method
!> We assume the initial guess of primary aqueous species is already set in the aqueous chemistry object
!> This subroutine is only used to read water types
subroutine initialise_conc_incr_coeff(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
    use aqueous_chemistry_m, only : aqueous_chemistry_c, int_array_c, LU_lin_syst, inf_norm_vec_real
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: icon(:)
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:)
    real(kind=8), intent(in) :: ctot(:)
    integer(kind=4), intent(out) :: niter !> number of iterations Newton-Raphson
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    
    real(kind=8), allocatable :: c1(:),c2_old(:),c2_new(:),conc_comp(:),c2_ig(:),log_c2k(:),log_c2(:),dc2_dc1(:,:),out_prod(:,:)
    real(kind=8), allocatable :: res(:) !> residual in Newton-Raphson
    real(kind=8), allocatable :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: abs_tol_res(:) !> absolute tolerances residues Newton-Raphson
    real(kind=8), allocatable :: mat_lin_syst(:,:),Se_aq_comp(:,:),K(:),u_aq(:),z2(:)
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:)
    logical :: flag_gas,flag_min,flag_wat,zero_flag
    type(int_array_c) :: indices_icon
    
    !external :: compute_res_init
!> Pre-process
    CV_flag=.false.
    if (size(icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) error stop
    
    if (sum(n_icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) error stop
    
    allocate(c2_old(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),&
        c2_new(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),&
        c2_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    ind_cstr=0 !> index of phase constrain
    call indices_icon%allocate_array(4)
    do i=1,indices_icon%num_cols
        call indices_icon%cols(i)%allocate_vector(n_icon(i))
    end do
    allocate(counters(4)) !> dim=n�icon options-2
    counters=0 !> counter for each modified icon option
    allocate(abs_tol_res(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        abs_tol_res(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol
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
            abs_tol_res(i)=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol
        else
            error stop "icon option not implemented"
        end if
    end do
!> Newton-Raphson
    allocate(res(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
    Jac_res(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(Delta_c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
    c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> chapuza

    c2_ig=1d-16 !> initial guess c2aq (chapuza)
    call this%set_conc_sec_aq_species(c2_ig) !> initial guess c2aq
    
    c1=1d0
    c1(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)=this%concentrations(1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) !> chapuza
    c2_old=1d-16
    call this%set_conc_sec_aq_species(c2_old(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)) !> initial guess c2aq
    !c2_old(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)=this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species) !> chapuza
    !c2_old(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)=ctot(indices_icon%cols(4)%col_1) !> chapuza
    niter=0
    this%log_act_coeffs=0d0 !> chapuza
!!> Squared charges species 
!    allocate(z2(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
!    z2=0d0 !> chapuza
!    z2(1:this%aq_phase%num_species)=this%aq_phase%z2
    do
    !> We update number of iterations
        niter=niter+1
        if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            error stop "Too many Newton iterations in initialisation"
        end if
        call this%compute_c2_from_c1_aq_Picard(c2_old,c2_new,niter_Picard,CV_flag)
        !call this%compute_ionic_act() !> we compute ionic activity
        !call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log activity coefficients aqueous species
        !call this%compute_activities()
        !call this%compute_log_act_coeff_wat()
        
        !print *, 10**this%log_act_coeffs
        !print *, c2_new
    !> Chapuza
        !c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)=this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)
    !> We compute resideu and Jacobian of residue using incremental coefficients
        call this%compute_res_Jac_res_incr_coeff(c2_new,indices_icon,n_icon,indices_constrains,ctot,res,Jac_res)
    !> We check convergence
        if (inf_norm_vec_real(res)<inf_norm_vec_real(abs_tol_res)) then !> CV reached
            CV_flag=.true.
            exit
        end if
    !> We solve linear system primary concentartions
        call LU_lin_syst(Jac_res,-res,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1)
        !> c1^(i+1)=c1^i+Delta_c1^i
        if (inf_norm_vec_real(Delta_c1/this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))<&
        this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            print *, inf_norm_vec_real(Delta_c1/this%concentrations(1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            error stop "Newton method not accurate enough in initialisation"
        end if
    !> We update primary concentrations
        call this%update_conc_aq_prim_species(Delta_c1,zero_flag)
    end do
    call this%compute_pH()
    call this%compute_salinity()
    call this%compute_alkalinity()

end subroutine