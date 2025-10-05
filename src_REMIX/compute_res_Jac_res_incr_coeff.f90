!> Computes residual and its  Jacobian using incremental coefficients
subroutine compute_res_Jac_res_incr_coeff(this,c2,indices_icon,n_icon,indices_constrains,ctot,res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, int_array_c
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c2(:) !> chapuza (dim=n_eq) (secondary concentrations for unperturbed residual)
    class(int_array_c), intent(in) :: indices_icon
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:)
    real(kind=8), intent(in) :: ctot(:)
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    
    real(kind=8), allocatable :: c1_pert(:),res_pert(:),c1(:),c2_pert(:),log_c2k(:),log_c2(:),dc2_dc1(:,:),out_prod(:,:)
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: abs_tol_res(:) !> absolute tolerances residues Newton-Raphson
    real(kind=8), allocatable :: mat_lin_syst(:,:),Se_aq_comp(:,:),K(:),u_aq(:),z2(:)
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr,niter
    integer(kind=4), allocatable :: indices_Jac(:),cols(:),ind_aq_species(:),counters(:)
    logical :: flag_gas,flag_min,flag_wat,CV_flag
    
    allocate(c1_pert(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),res_pert(&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    
    indices_Jac=indices_icon%get_vector_int() !> gets vector of icon indices
    call this%compute_res_init(indices_icon,n_icon,indices_constrains,ctot,res)
    c1=this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> unperturbed primary concentrations
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        !> we compute perturbation
        c1_pert=c1
        c1_pert(j)=c1_pert(j)+this%solid_chemistry%reactive_zone%CV_params%eps
        !> we set perturbed c1
        call this%set_conc_prim_species(c1_pert)
        !> we compute perturbed c2
        call this%compute_c2_from_c1_aq_Picard(c2,c2_pert,niter,CV_flag)
        !call this%compute_activities()
        !call this%compute_pH()
        !call this%compute_salinity()
        !call this%compute_alkalinity()
    !> We compute perturbed residual
        call this%compute_res_init(indices_icon,n_icon,indices_constrains,ctot,res_pert)
        Jac_res(:,j)=(res_pert-res)/this%solid_chemistry%reactive_zone%CV_params%eps
    end do
    call this%set_conc_prim_species(c1) !> unperturbed c1
    call this%set_conc_sec_aq_species(c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)) !> unperturbed c2aq
    call this%compute_ionic_act() !> we compute ionic activity
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log activity coefficients aqueous species
    call this%compute_activities_aq() !> chapuza
    call this%compute_log_act_coeff_wat()
end subroutine