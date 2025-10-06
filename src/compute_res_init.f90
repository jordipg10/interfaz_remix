!> Computes residual in initialisation
subroutine compute_res_init(this,indices_icon,n_icon,indices_constrains,ctot,res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst, int_array_c

    implicit none
    !> Arguments
    class(aqueous_chemistry_c) :: this
    class(int_array_c), intent(in) :: indices_icon
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:)
    real(kind=8), intent(in) :: ctot(:)
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
    
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr
    
    ind_cstr=0
    !> Concentrations
    do i=1,n_icon(1)
        res(indices_icon%cols(1)%col_1(i))=this%concentrations(indices_icon%cols(1)%col_1(i))-ctot(indices_icon%cols(1)%col_1(i))
    end do
    !> Aqueous components
    do i=1,n_icon(2)
        res(indices_icon%cols(2)%col_1(i))=dot_product(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(&
        indices_icon%cols(2)%col_1(i),1:this%aq_phase%num_species),this%concentrations)-ctot(indices_icon%cols(2)%col_1(i))
    end do
    !> Activities
    do i=1,n_icon(3)
        res(indices_icon%cols(3)%col_1(i))=(10**(this%log_act_coeffs(indices_icon%cols(3)%col_1(i))))*this%concentrations(&
        indices_icon%cols(3)%col_1(i))-ctot(indices_icon%cols(3)%col_1(i))
    end do
    !> Phases
    do i=1,n_icon(4)
        ind_cstr=ind_cstr+1 !> index constrain
        res(indices_icon%cols(4)%col_1(i))=dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(&
        indices_constrains(ind_cstr),1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species),&
        log10(this%activities(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)))+&
        log10(ctot(indices_icon%cols(4)%col_1(i)))-log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(&
        indices_constrains(ind_cstr))%eq_cst)
    end do
end subroutine