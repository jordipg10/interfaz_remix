!> Computes residual and its Jacobian analytically in initialisation (with cation exchange).
!> The residual vector and Jacobian are built in the water-type (quim_loc.dat) species ordering:
!> row/column k refer to the k-th species as listed in the water type definition. The chemical-system
!> matrices (Se, comp_mat_cst_act) are in chemical-system ordering, so their primary-species dimension
!> is mapped to the water ordering via this%ind_prim_species (ip(j) = water-type position of the j-th
!> chemical-system primary species). The last primary (position num_prim_species) is the exchange free site.
subroutine compute_res_Jac_res_anal_exch(this,conc,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,log_Jacobian_act_coeffs,CEC,&
    res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    use arrays_m, only: int_array_c
    use array_ops_m, only: get_ind_int !> map a chemical-system eq. reaction index to its per-water logK position
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: conc(:)
    class(int_array_c), intent(in) :: indices_icon
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> (:,1)=species position in water type, (:,2)=eq. reaction index in chem. system
    real(kind=8), intent(in) :: ctot(:)
    real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary (columns in chemical-system primary ordering)
    real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
    real(kind=8), intent(in) :: CEC
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson (water-type ordering)
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson (water-type ordering)

    real(kind=8), allocatable :: U_aq(:,:)
    integer(kind=4) :: i,j,ind_cstr,rxn,cp
    integer(kind=4) :: np,naqp,nsp,nsaq,naq
    real(kind=8) :: ln10

    ln10=log(1d1)
    ind_cstr=0 !> indices constrains
    Jac_res=0d0 !> Jacobian of residue
    associate (spc => this%solid_chemistry%reactive_zone%speciation_alg, &
               Se  => this%solid_chemistry%reactive_zone%chem_syst%Se, &
               U   => this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act, &
               ip  => this%ind_prim_species, &
               ier => this%solid_chemistry%reactive_zone%ind_eq_reacts)
    np  =spc%num_prim_species        !> primary species (aqueous primaries + cation-exchange free site)
    naqp=spc%num_aq_prim_species     !> aqueous primary species
    nsp =spc%num_species             !> all species
    nsaq=spc%num_sec_aq_species      !> secondary aqueous species
    naq =this%aq_phase%num_species   !> aqueous species (primary + secondary)
    !> Component matrix restricted to aqueous species (chemical-system column ordering)
    allocate(U_aq(naqp,naq))
    U_aq(:,1:naqp)=U(1:naqp,1:naqp)
    U_aq(:,naqp+1:naq)=U(1:naqp,np+1:naq)
!> Concentrations (icon=1): res_k = c_k - c_tot,k
    do i=1,n_icon(1)
        cp=indices_icon%cols(1)%col_1(i)
        res(cp)=conc(cp)-ctot(cp)
        Jac_res(cp,cp)=1d0
    end do
!> Aqueous components (icon=2): res_k = U_k . c_aq - u_tot,k
    do i=1,n_icon(2)
        cp=indices_icon%cols(2)%col_1(i)
        res(cp)=dot_product(U_aq(cp,:),this%concentrations([ip,this%ind_sec_species]))-ctot(cp)
        do j=1,naqp
            Jac_res(cp,ip(j))=U_aq(cp,j)+dot_product(U_aq(cp,naqp+1:naq),dc2_dc1(1:nsaq,j))
        end do
    end do
!> Activities (icon=3): res_k = a_k - a_tot,k
    do i=1,n_icon(3)
        cp=indices_icon%cols(3)%col_1(i)
        res(cp)=this%activities(cp)-ctot(cp)
        Jac_res(cp,cp)=1d1**this%log_act_coeffs(cp)
    end do
!> Phases (icon=4): res_k = Se . log10(a_prim) + log10(ctot) - logK
    do i=1,n_icon(4)
        ind_cstr=ind_cstr+1 !> index constrain
        cp=indices_icon%cols(4)%col_1(i)
        rxn=indices_constrains(ind_cstr,2) !> eq. reaction index in chemical system
        res(cp)=dot_product(Se(rxn,1:naqp),log10(this%activities(ip)))+&
            log10(ctot(cp))-this%logK(get_ind_int(rxn,ier))
        do j=1,naqp
            Jac_res(cp,ip(j))=Se(rxn,j)/(this%concentrations(ip(j))*ln10)
        end do
    end do
!> Cation exchange capacity (free-site primary, position np): res_np = U_np . c - CEC
!> The free-site column (np) couples the exchange secondaries via dc2_dc1(:,np).
    res(np)=dot_product(U(np,1:naqp),conc(ip))+dot_product(U(np,np:nsp),conc(np:nsp))-CEC
    do j=1,naqp
        Jac_res(np,ip(j))=U(np,j)+dot_product(U(np,np+1:nsp),dc2_dc1(:,j))
    end do
    Jac_res(np,np)=U(np,np)+dot_product(U(np,np+1:nsp),dc2_dc1(:,np))
    end associate
end subroutine