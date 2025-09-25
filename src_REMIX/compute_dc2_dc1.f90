!> Computes Jacobian of secondary concentrations with respect to primary concentrations
!! \f$\partial c_2 / \partial c_1\f$
subroutine compute_dc2_dc1(this,out_prod,c1,c2,dc2_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, diag_matrix_c, id_matrix

    implicit none
!> Variables
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_d_I and z^2
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(in) :: c2(:) !> chapuza (dim=n_eq)
    real(kind=8), intent(out) :: dc2_dc1(:,:) !> Jacobian of c2 with respect to c1
    
    integer(kind=4) :: i,j
    type(diag_matrix_c) :: diag_c2,diag_c1
    real(kind=8), allocatable :: mat_lin_syst(:,:),aux_mat(:,:),indep_mat(:,:),indep_term_aux(:),log_dc2_dc1(:,:)
    
    call diag_c2%set_diag_matrix(c2)
    call diag_c1%set_diag_matrix(c1)
!> Then we compute log-Jacobian activity coefficients with respect to secondary aqueous concentrations
    !call this%aq_phase%compute_d_log_gamma_d_c2_aq(this%d_log_gamma_d_I,this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species),this%log_Jacobian_act_coeffs(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species))
    !d_log_gamma_d_log_c=0d0
    !d_log_gamma_d_log_c(1:this%aq_phase%num_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)=this%log_Jacobian_act_coeffs(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)
    !d_log_gamma1_d_log_c2=d_log_gamma_d_log_c(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species)
    !d_log_gamma2_d_log_c2=d_log_gamma_d_log_c(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species)
!> We compute linear system with logarithms
    allocate(log_dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    aux_mat=5d-1*log(1d1)*(matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,out_prod(1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_species))-&
    out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_species))
    !> Finally, we obtain the matrix that premultiplies the log-Jacobian
    mat_lin_syst=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-diag_c2%prod_mat_diag_mat(aux_mat)
    !> We compute independent matrix
    aux_mat=5d-1*log(1d1)*diag_c1%prod_mat_diag_mat(out_prod(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
    1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))+&
    id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    indep_mat=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,aux_mat)-5d-1*log(1d1)*&
    diag_c1%prod_mat_diag_mat(out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_species,1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    !> For each primary species, we solve linear system
    !allocate(indep_mat(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !!> Independent term j-th primary species
    !    indep_term_aux=this%aq_phase%z2(j)*5d-1*log(1d1)*this%concentrations(j)*this%d_log_gamma_d_I(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    !    indep_term_aux(j)=indep_term_aux(j)+1d0
    !    indep_term(:,j)=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,indep_term_aux)-this%aq_phase%z2(j)*5d-1*log(1d1)*this%concentrations(j)*this%d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species)
     !> Linear system solver
        call LU_lin_syst(mat_lin_syst,indep_mat(:,j),this%solid_chemistry%reactive_zone%CV_params%zero,log_dc2_dc1(:,j))
    end do
    !> Finally, we compute dc2_dc1
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2_dc1(i,j)=c2(i)*log_dc2_dc1(i,j)/c1(j)
        end do
    end do
end subroutine 