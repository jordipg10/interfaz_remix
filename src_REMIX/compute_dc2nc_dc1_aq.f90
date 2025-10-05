!> Computes Jacobian of secondary variable activity species concentrations with respect to primary species
!! \f$\partial c_{2,nc} / \partial c_1\f$
!! We assume all primary species are aqueous
subroutine compute_dc2nc_dc1_aq(this,c2nc,out_prod,dc2nc_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, diag_matrix_c, id_matrix
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c2nc(:)
    real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_nc_d_I and z_nc^2
    real(kind=8), intent(out) :: dc2nc_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j
    type(diag_matrix_c) :: diag_c2nc,diag_c1
    real(kind=8), allocatable :: aux_mat(:,:),mat_lin_syst(:,:),indep_mat(:,:)
    real(kind=8), allocatable :: log_dc2nc_dc1(:,:),d_log_gamma2nc_d_log_c2nc(:,:)
    
!> We compute linear system with logarithms
    allocate(log_dc2nc_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    
    call diag_c2nc%set_diag_matrix(c2nc)
    aux_mat=5d-1*log(1d1)*(matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,out_prod(1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species))-out_prod(&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    !> We obtain the matrix that premultiplies the log-Jacobian
    mat_lin_syst=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-diag_c2nc%prod_mat_diag_mat(aux_mat)
    !> We compute independent matrix
    call diag_c1%set_diag_matrix(this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    aux_mat=5d-1*log(1d1)*diag_c1%prod_mat_diag_mat(out_prod(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
    1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))+&
    id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    indep_mat=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,aux_mat)-5d-1*log(1d1)*&
    diag_c1%prod_mat_diag_mat(out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) 
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> Linear system solver
        call LU_lin_syst(mat_lin_syst,indep_mat(:,j),this%solid_chemistry%reactive_zone%CV_params%zero,log_dc2nc_dc1(:,j))
    end do
    !> Finally, we compute dc2nc_dc1
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2nc_dc1(i,j)=c2nc(i)*log_dc2nc_dc1(i,j)/this%concentrations(j)
        end do
    end do 
end subroutine 