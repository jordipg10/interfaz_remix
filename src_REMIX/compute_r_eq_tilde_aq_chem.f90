!!> Computes equilibrium reaction rates from secondary aqueous variable activity species concentrations
!subroutine compute_r_eq_tilde_aq_chem(this,this_old,mixing_ratios,mixing_waters,Delta_t,porosity,r_eq_tilde)
!    use aqueous_chemistry_m
!    use metodos_sist_lin_m
!    implicit none
!    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at time step k+1
!    !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
!    class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at time step k (nombre muy malo)
!    real(kind=8), intent(in) :: mixing_ratios(:)
!    class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
!    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
!    real(kind=8), intent(in) :: porosity
!    real(kind=8), intent(out) :: r_eq_tilde(:)
!    
!    integer(kind=4) :: i,n,niter,rk_ind,l
!    real(kind=8), allocatable :: c2nc_tilde(:),A(:,:),b(:),c2nc(:),R_eq(:)
!    
!    allocate(c2nc_tilde(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species))
!!> We compute c2nc tilde
!    call this_old%compute_c2nc_tilde_aq_chem(mixing_ratios,mixing_waters,c2nc_tilde)
!!> Mínimos cuadrados a partir de balance masas
!    A=matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),transpose(this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)))
!    c2nc=this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)
!    b=matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),c2nc-c2nc_tilde)
!    allocate(R_eq(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)) !> R_eq=Delta_t*r_eq*Lambda^T_j/phi
!    if (inf_norm_vec_real(b)<this%solid_chemistry%reactive_zone%CV_params%zero) then
!        R_eq=0d0
!    else
!        call LU_lin_syst(A,b,this%solid_chemistry%reactive_zone%CV_params%zero,R_eq) !> linear system solver
!    end if
!    r_eq_tilde=R_eq*porosity/Delta_t !> r_eq_tilde_j=r_eq*Lambda^T_j
!end subroutine