!> Computes residual and its Jacobian analytically in initialisation
subroutine compute_res_Jac_res_anal_exch(this,conc,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,log_Jacobian_act_coeffs,CEC,&
    res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, int_array_c
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: conc(:)
    class(int_array_c), intent(in) :: indices_icon
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:)
    real(kind=8), intent(in) :: ctot(:)
    real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
    real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
    real(kind=8), intent(in) :: CEC
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    
    real(kind=8), allocatable :: c1(:),c2(:),c2nc_init(:),c1k(:),c2k(:),log_c2k(:),log_c2(:),out_prod(:,:)
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: abs_tol_res(:) !> absolute tolerances residues Newton-Raphson
    real(kind=8), allocatable :: mat_lin_syst(:,:),U_aq(:,:),K(:),z2(:)
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:)
    logical :: flag_gas,flag_min,flag_wat
    
    ind_cstr=0 !> indices constrains
    Jac_res=0d0 !> Jacobian of residue
    allocate(U_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,this%aq_phase%num_species))
    U_aq(:,1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)= &
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species, &
        1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
    U_aq(:,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)= &
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species, &
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)
!> Concentrations
    do i=1,n_icon(1)
        res(indices_icon%cols(1)%col_1(i))=conc(indices_icon%cols(1)%col_1(i))-ctot(indices_icon%cols(1)%col_1(i))
        Jac_res(indices_icon%cols(1)%col_1(i),indices_icon%cols(1)%col_1(i))=1d0
    end do
!> Aqueous components
    do i=1,n_icon(2)
        res(indices_icon%cols(2)%col_1(i))=dot_product(U_aq(indices_icon%cols(2)%col_1(i),:),conc)-&
        ctot(indices_icon%cols(2)%col_1(i))
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            Jac_res(indices_icon%cols(2)%col_1(i),j)=U_aq(indices_icon%cols(2)%col_1(i),j)+dot_product(U_aq(&
            indices_icon%cols(2)%col_1(i), &
            this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species),&
            dc2_dc1(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species,j))
        end do
    end do
!> Activities
    do i=1,n_icon(3)
        res(indices_icon%cols(3)%col_1(i))=(10**(this%log_act_coeffs(indices_icon%cols(3)%col_1(i))))*&
        conc(indices_icon%cols(3)%col_1(i))-ctot(indices_icon%cols(3)%col_1(i))
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            Jac_res(indices_icon%cols(3)%col_1(i),j)=(10**(this%log_act_coeffs(indices_icon%cols(3)%col_1(i))))*(&
            log_Jacobian_act_coeffs(indices_icon%cols(3)%col_1(i),j)+dot_product(log_Jacobian_act_coeffs(&
            indices_icon%cols(3)%col_1(i),this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species),dc2_dc1(:,j)))*(conc(indices_icon%cols(3)%col_1(i))/&
            conc(j))
        end do
        Jac_res(indices_icon%cols(3)%col_1(i),indices_icon%cols(3)%col_1(i))=Jac_res(indices_icon%cols(3)%col_1(i),&
        indices_icon%cols(3)%col_1(i)) + 10**(this%log_act_coeffs(indices_icon%cols(3)%col_1(i)))
    end do
!> Phases
    do i=1,n_icon(4)
        ind_cstr=ind_cstr+1 !> index constrain
        res(indices_icon%cols(4)%col_1(i))=dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(&
        indices_constrains(ind_cstr),1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
        log10(this%activities(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)))+&
        log10(ctot(indices_icon%cols(4)%col_1(i)))-log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(&
        indices_constrains(ind_cstr))%eq_cst)
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            Jac_res(indices_icon%cols(4)%col_1(i),j)=dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(&
            indices_constrains(ind_cstr),1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),&
            log_Jacobian_act_coeffs(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,j)+&
            matmul(log_Jacobian_act_coeffs(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species),dc2_dc1(:,j)))
            Jac_res(indices_icon%cols(4)%col_1(i),j)=Jac_res(indices_icon%cols(4)%col_1(i),j)+&
            this%solid_chemistry%reactive_zone%chem_syst%Se(indices_constrains(ind_cstr),j)
            Jac_res(indices_icon%cols(4)%col_1(i),j)=Jac_res(indices_icon%cols(4)%col_1(i),j)/(conc(j)*log(1d1))
        end do
    end do
!> Cation exchange capacity
    res(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=dot_product(&
    this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,:),conc)-CEC
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        Jac_res(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,j)=&
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,j)+dot_product(&
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_species),dc2_dc1(:,j))
    end do
end subroutine