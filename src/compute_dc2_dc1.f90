!> \file compute_dc2_dc1.f90
!> \brief Computes Jacobian of secondary species concentrations with respect to primary species concentrations including activity coefficient effects
!> \details This subroutine calculates the Jacobian matrix ∂c₂/∂c₁ for non-ideal solutions where activity 
!> coefficients depend on ionic strength. The computation involves solving a linear system in logarithmic 
!> variables and then transforming back to linear concentrations.
!>
!> The Jacobian relates changes in secondary species concentrations to changes in primary species 
!> concentrations, accounting for:
!> - Mass action law relationships via stoichiometric matrix Sₑ₁*
!> - Activity coefficient dependencies on ionic strength
!> - Ionic strength dependencies on all species concentrations
!>
!> Mathematical formulation in logarithmic form:
!> \f[
!> \frac{\partial \log c_2}{\partial \log c_1} = \text{solution of linear system}
!> \f]
!> Then convert to linear form:
!> \f[
!> \frac{\partial c_2}{\partial c_1} = \frac{c_2}{c_1} \cdot \frac{\partial \log c_2}{\partial \log c_1}
!> \f]
!>
!> The linear system solved is:
!> \f[
!> \left[\mathbf{I} - \mathbf{D}_{c_2} \cdot \mathbf{A}_{act}\right] \cdot \frac{\partial \log c_2}{\partial \log c_1} = \mathbf{B}_{indep}
!> \f]
!> where activity coefficient corrections are incorporated via the outer product of ∂(log γ)/∂I and z².
!>
!> \param[in,out] this Aqueous chemistry object containing reactive zone, speciation algebra, and species data
!> \param[in] out_prod Outer product matrix: [∂(log γ)/∂I] ⊗ [z²] for activity coefficient corrections [-]
!> \param[in] c1 Primary species concentrations [C]
!> \param[in] c2 Secondary species concentrations [C]
!> \param[out] dc2_dc1 Jacobian matrix: ∂c₂/∂c₁ relating secondary to primary concentrations [-]

subroutine compute_dc2_dc1(this,out_prod,c1,c2,dc2_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use arrays_m, only: id_matrix, diag_matrix_c
    use metodos_sist_lin_m, only: LU_lin_syst
    implicit none !> Enforce explicit variable declarations
!> Variables
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing reactive zone, aqueous phase, speciation algebra, and species concentrations [-]
    real(kind=8), intent(in) :: out_prod(:,:) !> Outer product matrix: [∂(log γ)/∂I] ⊗ [z²] for activity coefficient corrections in Jacobian [-]
    real(kind=8), intent(in) :: c1(:) !> Primary species concentrations (dimension = n_p) [C]
    real(kind=8), intent(in) :: c2(:) !> Secondary species concentrations (dimension = n_eq) [C]
    real(kind=8), intent(out) :: dc2_dc1(:,:) !> Jacobian matrix: ∂c₂/∂c₁ relating changes in secondary to changes in primary concentrations [-]
    
    integer(kind=4) :: i !> Loop counter for secondary species (rows in Jacobian) [-]
    integer(kind=4) :: j !> Loop counter for primary species (columns in Jacobian) [-]
    type(diag_matrix_c) :: diag_c2 !> Diagonal matrix with secondary species concentrations on diagonal [-]
    type(diag_matrix_c) :: diag_c1 !> Diagonal matrix with primary species concentrations on diagonal [-]
    real(kind=8), allocatable :: mat_lin_syst(:,:) !> Linear system matrix: I - D_c2·A_act for solving log-Jacobian [-]
    real(kind=8), allocatable :: aux_mat(:,:) !> Auxiliary matrix for intermediate calculations in activity coefficient corrections [-]
    real(kind=8), allocatable :: indep_mat(:,:) !> Independent term matrix (right-hand side) for linear system: one column per primary species [-]
    real(kind=8), allocatable :: indep_term_aux(:) !> Auxiliary vector for independent term (not used in current implementation) [-]
    real(kind=8), allocatable :: log_dc2_dc1(:,:) !> Log-Jacobian matrix: ∂(log c₂)/∂(log c₁) in logarithmic variables [-]
    
    call diag_c2%set_diag_matrix(c2) !> Create diagonal matrix with secondary concentrations: D_c2 = diag(c₂)
    call diag_c1%set_diag_matrix(c1) !> Create diagonal matrix with primary concentrations: D_c1 = diag(c₁)
!> Then we compute log-Jacobian activity coefficients with respect to secondary aqueous concentrations (COMMENTED OUT - alternative approach)
    !call this%aq_phase%compute_d_log_gamma_d_c2_aq(this%d_log_gamma_d_I,this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species),this%log_Jacobian_act_coeffs(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)) !> Commented: Direct computation of ∂(log γ)/∂c₂
    !d_log_gamma_d_log_c=0d0 !> Commented: Initialize full log-Jacobian of activity coefficients
    !d_log_gamma_d_log_c(1:this%aq_phase%num_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species)=this%log_Jacobian_act_coeffs(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species) !> Commented: Assign secondary species block
    !d_log_gamma1_d_log_c2=d_log_gamma_d_log_c(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species) !> Commented: Extract ∂(log γ₁)/∂(log c₂) block
    !d_log_gamma2_d_log_c2=d_log_gamma_d_log_c(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species) !> Commented: Extract ∂(log γ₂)/∂(log c₂) block
!> We compute linear system matrix with logarithmic variables and activity coefficient corrections
    allocate(log_dc2_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,& !> Allocate log-Jacobian matrix with dimensions: n_eq × n_p
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> Number of columns = number of primary species
    aux_mat=5d-1*log(1d1)*(matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,out_prod(1:& !> Compute activity coefficient correction matrix: (ln(10)/2)·[Sₑ₁*·OP₁₂ - OP₂₂], start with Sₑ₁* times outer product block (primary→secondary)
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,& !> Rows: primary species indices in outer product
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Columns start: first secondary species
    this%solid_chemistry%reactive_zone%speciation_alg%num_species))-& !> Columns end: last species, then subtract next term
    out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Subtract outer product block (secondary→secondary): rows start at first secondary
    this%solid_chemistry%reactive_zone%speciation_alg%num_species,& !> Rows end: last species
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Columns start: first secondary species
    this%solid_chemistry%reactive_zone%speciation_alg%num_species)) !> Columns end: last species; multiply by ln(10)/2 for log₁₀ to ln conversion factor
    !> Finally, we obtain the matrix that premultiplies the log-Jacobian in linear system
    mat_lin_syst=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-diag_c2%prod_mat_diag_mat(aux_mat) !> Form left-hand side matrix: I - D_c2·A_act where A_act includes activity corrections
    !> We compute independent matrix (right-hand side) incorporating activity effects from primary species
    aux_mat=5d-1*log(1d1)*diag_c1%prod_mat_diag_mat(out_prod(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,& !> Compute (ln(10)/2)·D_c1·OP₁₁: activity corrections for primary-primary interactions
    1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))+& !> Outer product block: primary→primary, multiply by diagonal concentration matrix
    id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> Add identity matrix I for direct concentration dependencies
    indep_mat=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,aux_mat)-5d-1*log(1d1)*& !> Compute independent matrix: Sₑ₁*·[...] - (ln(10)/2)·D_c1·OP₂₁
    diag_c1%prod_mat_diag_mat(out_prod(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Subtract activity correction term: D_c1 times outer product block (secondary→primary)
    this%solid_chemistry%reactive_zone%speciation_alg%num_species,1:& !> Rows: secondary species in outer product
    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> Columns: primary species in outer product
    !> For each primary species, we solve linear system: mat_lin_syst·∂(log c₂)/∂(log c₁ⱼ) = indep_mat(:,j)
    !allocate(indep_mat(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)) !> Commented: Alternative allocation approach
    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Loop over all primary species (each column of Jacobian)
    !!> Independent term j-th primary species (COMMENTED OUT - alternative formulation)
    !    indep_term_aux=this%aq_phase%z2(j)*5d-1*log(1d1)*this%concentrations(j)*this%d_log_gamma_d_I(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> Commented: Compute activity contribution for j-th primary on primary species
    !    indep_term_aux(j)=indep_term_aux(j)+1d0 !> Commented: Add direct dependence (∂ log cⱼ/∂ log cⱼ = 1)
    !    indep_term(:,j)=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,indep_term_aux)-this%aq_phase%z2(j)*5d-1*log(1d1)*this%concentrations(j)*this%d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_species) !> Commented: Apply stoichiometry and subtract secondary contributions
     !> Linear system solver: Solve for ∂(log c₂)/∂(log c₁ⱼ) for j-th primary species
        call LU_lin_syst(mat_lin_syst,indep_mat(:,j),this%solid_chemistry%reactive_zone%CV_params%zero,log_dc2_dc1(:,j)) !> Solve linear system using LU decomposition: A·x = b where x = log-Jacobian column j
    end do !> End loop over primary species
    !> Finally, we compute dc2_dc1 by converting from logarithmic to linear Jacobian
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !> Loop over all secondary species (rows in Jacobian)
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Loop over all primary species (columns in Jacobian)
            dc2_dc1(i,j)=c2(i)*log_dc2_dc1(i,j)/c1(j) !> Convert log-Jacobian to linear Jacobian: ∂c₂ᵢ/∂c₁ⱼ = (c₂ᵢ/c₁ⱼ)·∂(log c₂ᵢ)/∂(log c₁ⱼ)
        end do !> End loop over primary species
    end do !> End loop over secondary species
end subroutine 