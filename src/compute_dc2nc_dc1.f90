!> Computes Jacobian of secondary variable activity species concentrations
!> \f$\partial c_{2,nc} / \partial c_1\f$
subroutine compute_dc2v_dc1(this,c1,c2v,dc2v_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real, outer_prod_vec
    use metodos_sist_lin_m, only: LU_lin_syst
    use arrays_m, only: diag_matrix_c, id_matrix
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza
    real(kind=8), intent(in) :: c2v(:) !> chapuza
    !real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_nc_d_I and z_nc^2
    real(kind=8), intent(out) :: dc2v_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j,n_v
    type(diag_matrix_c) :: diag_c1,diag_c2v
    real(kind=8), allocatable :: aux_mat(:,:),mat_lin_syst(:,:),indep_mat(:,:),log_Jacobian(:,:),log_dc2v_dc1(:,:),d_log_gamma_d_I(:),out_prod(:,:)
    
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(log_dc2v_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(mat_lin_syst(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    allocate(indep_mat(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(aux_mat(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    !allocate(d_log_gamma_d_I(n_v))
    !allocate(out_prod(n_v,n_v))
    
    !d_log_gamma_d_I=0d0 !> Initialize activity coefficient derivatives to zero (temporary workaround)
    d_log_gamma_d_I=this%compute_d_log_gamma_d_I_aq_chem() !> Compute ∂(log γᵢ)/∂I for all variable activity species (from Debye-Hückel or extended models)

    out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:n_v)) !> Compute outer product: [∂(log γ)/∂I] ⊗ [z²] for activity corrections in Jacobian
    !> We compute linear system with logarithms
    
    
    call diag_c2v%set_diag_matrix(c2v)
    aux_mat=5d-1*log(1d1)*(matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,out_prod(1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))-out_prod(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    !> We obtain the matrix that premultiplies the log-Jacobian
    mat_lin_syst=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-diag_c2v%prod_mat_diag_mat(aux_mat)
    !> We compute independent matrix
    call diag_c1%set_diag_matrix(c1)
    aux_mat=5d-1*log(1d1)*diag_c1%prod_mat_diag_mat(out_prod(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
        1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))+&
        id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    indep_mat=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,aux_mat)-5d-1*log(1d1)*&
        diag_c1%prod_mat_diag_mat(out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) 
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> Linear system solver
        call LU_lin_syst(mat_lin_syst,indep_mat(:,j),this%solid_chemistry%reactive_zone%CV_params%zero,log_dc2v_dc1(:,j))
    end do
    !> Finally, we compute dc2v_dc1
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2v_dc1(i,j)=c2v(i)*log_dc2v_dc1(i,j)/c1(j)
        end do
    end do
end subroutine