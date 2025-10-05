!> Computes Jacobian of Newton residual with respect to primary concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has equilibrium and kinetic reactions
!> We assume all primary species are aqueous
!> We apply lumping in this subroutine
subroutine compute_dfk_dc1_aq_EfI_lump(this,c2nc,drk_dc,Delta_t,dfk_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, outer_prod_vec
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c2nc(:)
    real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates
    !real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: dfk_dc1(:,:) !> Jacobian Newton residual - primary concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: dc2nc_dc1(:,:) !> Jacobian secondary variable activity concentrations - primary concentrations
    real(kind=8), allocatable :: drk_dc1(:,:) !> Jacobian kinetic reaction rates - primary concentrations
    real(kind=8), allocatable :: drk_dc2nc(:,:) !> Jacobian kinetic reaction rates - secondary variable activity concentrations
    real(kind=8), allocatable :: d_log_gamma_d_I(:) !> derivative of log_10 activity coefficients of variable activity species with respect to ionic activity
    real(kind=8), allocatable :: out_prod(:,:) !> Outer product between d_log_gamma_d_I and squared charges of variable activity species
!> Pre-process
    allocate(dc2nc_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
!> Process
    !> First we compute d_log_gamma_d_I
        call this%compute_d_log_gamma_d_I_aq_chem(d_log_gamma_d_I)
    !> Outer product
        out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    !> Jacobian of secondary variable activity concentrations
        call this%compute_dc2nc_dc1_aq(c2nc,out_prod,dc2nc_dc1) 
    !> We separate primary and secondary variable activity species in Jacobian of kinetic reaction rates       
        drk_dc1=drk_dc(:,1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
        drk_dc2nc=drk_dc(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)
    !> We compute Jacobian Newton residual
        dfk_dc1=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+matmul(&
            this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),dc2nc_dc1)-Delta_t*&
            matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc1+matmul(drk_dc2nc,dc2nc_dc1))
!> Post-process
    deallocate(d_log_gamma_d_I,dc2nc_dc1,drk_dc1,drk_dc2nc,out_prod)
end subroutine