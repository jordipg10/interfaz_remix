!> Computes residual and its Jacobian analytically in initialisation assuming ideal activities
subroutine compute_res_Jac_res_anal_ideal(this,indices_icon,n_icon,indices_constrains,ctot,dc2aq_dc1,res,Jac_res)
use aqueous_chemistry_m, only: aqueous_chemistry_c
use matrices_m, only: int_array_c
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    class(int_array_c), intent(in) :: indices_icon !> indices of icon options
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> indices of constrains in chemical system
    real(kind=8), intent(in) :: ctot(:) !< must be ordered as in input file
    real(kind=8), intent(in) :: dc2aq_dc1(:,:) !> Jacobian aqueous secondary-aqueous primary
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
    
    real(kind=8), allocatable :: c1(:),c2(:),c2nc_init(:),c1k(:),c2k(:),log_c2k(:),log_c2(:),out_prod(:,:)
    real(kind=8), allocatable :: Delta_c1(:) !> c1^(i+1)-c1^i (Newton)
    real(kind=8), allocatable :: abs_tol_res(:) !> absolute tolerances residues Newton-Raphson
    real(kind=8), allocatable :: mat_lin_syst(:,:),U_aq(:,:),K(:),z2(:)
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:)
    logical :: flag_gas,flag_min,flag_wat
    
    ind_cstr=0
    Jac_res=0d0 !> Jacobian residue
    allocate(U_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
        this%aq_phase%num_species))
    U_aq(:,1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)=&
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
    U_aq(:,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)=&
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%aq_phase%num_species)
    !> Concentrations
    !print *, indices_icon%cols(1)%col_1
    do i=1,n_icon(1)
        res(indices_icon%cols(1)%col_1(i))=this%concentrations(indices_icon%cols(1)%col_1(i))-ctot(indices_icon%cols(1)%col_1(i))
        Jac_res(indices_icon%cols(1)%col_1(i),indices_icon%cols(1)%col_1(i))=1d0
    end do
    !> Aqueous components
    !print *, this%ind_prim_species
    !print *, this%ind_sec_species
    do i=1,n_icon(2)
        res(indices_icon%cols(2)%col_1(i))=dot_product(U_aq(indices_icon%cols(2)%col_1(i),:),this%concentrations(&
            [this%ind_prim_species,this%ind_sec_species]))-& !> chapuza
            ctot(indices_icon%cols(2)%col_1(i))
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            Jac_res(indices_icon%cols(2)%col_1(i),this%ind_prim_species(j))=U_aq(indices_icon%cols(2)%col_1(i),j)+dot_product(U_aq(&
                indices_icon%cols(2)%col_1(i),this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
                this%aq_phase%num_species),dc2aq_dc1(:,j))
        end do
    end do
    !> Activities
    do i=1,n_icon(3)
        res(indices_icon%cols(3)%col_1(i))=this%activities(indices_icon%cols(3)%col_1(i))-ctot(indices_icon%cols(3)%col_1(i))
        Jac_res(indices_icon%cols(3)%col_1(i),indices_icon%cols(3)%col_1(i))=10**this%log_act_coeffs(indices_icon%cols(3)%col_1(i))
    end do
    !> Phases
    !print *, indices_icon%cols(4)%col_1
    do i=1,n_icon(4)
        ind_cstr=ind_cstr+1 !> index constrain
        res(indices_icon%cols(4)%col_1(i))=dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(&
            indices_constrains(ind_cstr,2),1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),&
            log10(this%activities(this%ind_prim_species)))+log10(ctot(indices_icon%cols(4)%col_1(i)))-&
            log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr,2))%eq_cst)
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            Jac_res(indices_icon%cols(4)%col_1(i),this%ind_prim_species(j))=this%solid_chemistry%reactive_zone%chem_syst%Se(&
                indices_constrains(ind_cstr,2),j)/(this%concentrations(this%ind_prim_species(j))*log(1d1))
        end do
    end do
end subroutine