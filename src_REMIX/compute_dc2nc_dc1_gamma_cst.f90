!> Computes Jacobian of secondary variable activity species concentrations assuming constant activity coefficients
!> \f$\partial c_{2,nc} / \partial c_1\f$ with \f$gamma\f$ constant
subroutine compute_dc2nc_dc1_gamma_cst(this,dc2nc_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, diag_matrix_c

    implicit none
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(out) :: dc2nc_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j,n_p,n_e,n_nc_aq
    type(diag_matrix_c) :: id_c1,id_c2
    real(kind=8), allocatable :: aux_mat(:,:),c1(:),c2nc(:)
    
!> Process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_nc_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species

    c1=this%get_c1()
    call id_c1%set_diag_matrix(1d0/c1)
    aux_mat=this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star
    do i=1,n_e
        do j=1,n_p
            aux_mat(i,j)=aux_mat(i,j)*id_c1%diag(j)
        end do
    end do
    c2nc=this%get_c2nc()
    call id_c2%set_diag_matrix(c2nc)
    !id_c2=diag(c2nc)
    dc2nc_dc1=aux_mat
    do i=1,n_e
        do j=1,n_p
            dc2nc_dc1(i,j)=dc2nc_dc1(i,j)*id_c2%diag(i)
        end do
    end do
end subroutine 