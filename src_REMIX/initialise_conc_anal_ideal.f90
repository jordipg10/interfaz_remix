!> Computes species concentrations with data from water type definition using Newton method
!> We assume the initial guess of primary aqueous species is already set in the aqueous chemistry object
!> We assume ideal conditions & no exchange reactions
!! We assume concentrations are given in molality
!> This subroutine is only used to read water types based in CHEPROO
subroutine initialise_conc_anal_ideal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, int_array_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: icon(:) !> initial condition type
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> indices of constrains in reactive zone stoichiometric matrix
    real(kind=8), intent(in) :: ctot(:) !> data given (must be ordered as in "aq_phase" attribute)
    integer(kind=4), intent(out) :: niter !> number of iterations Newton-Raphson
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    !character(len=*), intent(in), optional :: path_py !> path output Python
!> Variables
    real(kind=8), allocatable :: c1(:),c2(:),conc_comp(:),a1(:),log_gamma2(:),log_c2(:),dc2_dc1(:,:)
    real(kind=8), allocatable :: res(:) !> residual in Newton-Raphson
    real(kind=8), allocatable :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: tol_res(:) !> tolerance residues Newton-Raphson
    real(kind=8), allocatable ::mat_lin_syst(:,:),Se_aq_comp(:,:),K(:),u_aq(:),z2(:),d_log_gamma_d_I(:)
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:)
    logical :: CV_flag_Picard
    type(int_array_c) :: indices_icon
    
!> Pre-process
    CV_flag=.false.
    if (size(icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
        error stop "Dimension error in icon"
    else if (sum(n_icon)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
         error stop "Dimension error in n_icon"
    end if
    allocate(c2(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
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
            indices_icon%cols(4)%col_1(counters(4))=indices_constrains(counters(4),1)
            tol_res(i)=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol !> logarithmic absolute tolerance
        else
            error stop "icon option not implemented"
        end if
    end do
!> Newton-Raphson
    allocate(res(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),&
        Jac_res(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    allocate(Delta_c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
    allocate(dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))    
    !allocate(c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    
    niter=0
    call this%set_act_aq_species() !> we set ideal activities of aqueous species
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    !call this%compute_salinity()
    !call this%compute_molarities()
    log_gamma2=this%get_log_gamma2() !> we get log_10 activity coefficients of secondary species
    do
        !> We update number of iterations
        niter=niter+1
        if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            error stop "Too many Newton iterations in initialisation"
        end if
        call this%compute_activities_diss_solids() !> chapuza
        a1=this%get_a1() !> chapuza
        call this%compute_c2_from_a1_ideal(a1,log_gamma2,c2)
        c1=this%get_c1_aq() !> we get primary aqueous concentrations (chapuza)
        !> We compute Jacobian secondary-primary concentrations
        call this%compute_dc2_dc1_ideal(c1,c2,dc2_dc1)
        !> We compute residual and its Jacobian
        call this%compute_res_Jac_res_anal_ideal(indices_icon,n_icon,indices_constrains,ctot,dc2_dc1(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species,:),res,Jac_res)
        !> We check convergence
        if (inf_norm_vec_real(res)<inf_norm_vec_real(tol_res)) then !> CV reached
            CV_flag=.true.
            exit
        end if
        !> We solve linear system aqueous primary concentartions
        call LU_lin_syst(Jac_res,-res,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1)
        !> c1^(i+1)=c1^i+Delta_c1^i
        if (inf_norm_vec_real(Delta_c1/c1)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            print *, inf_norm_vec_real(Delta_c1/c1)
            error stop "Newton method not accurate enough in initialisation"
        end if
        !> We update primary concentrations
        call this%update_conc_aq_prim_species_bis(Delta_c1(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    end do
    call this%compute_pH()
    call this%compute_activities_diss_solids()
    call this%compute_salinity() !> we compute salinity to change units
    !call this%compute_molarities() !> we change units to molarities
    call this%compute_alkalinity()
    !call this%compute_activities_diss_solids()
!> Post-process
    if (associated(this%gas_chemistry)) then
        nullify(this%gas_chemistry)
    end if
end subroutine