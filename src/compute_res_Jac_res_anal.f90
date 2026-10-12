!> Computes residual and its Jacobian analytically in initialisation (non-ideal water types).
!!
!! Unknowns of the Newton method: the aqueous primary concentrations c1 in STOICHIOMETRIC order (the
!! order of Se, Se_1_star and comp_mat_cst_act, i.e. of get_c1 and of the columns of dc2aq_dc1); primary
!! j is stored at concentration index p(j) = ind_prim_species(j), aqueous secondary i at
!! s(i) = ind_sec_species(i). The constraints (icon, ctot, indices_icon) are given by INPUT position,
!! which is the concentration index of the constrained primary (read_wat_type_CHEPROO stores
!! concentrations(k) = guess, icons(k), ctots(k)); the row of constraint k is the primary jk with p(jk) = k,
!! so that res and Jac_res come out in the unknowns' order and the Newton update Delta_c1 is applied to the
!! right species (the previous version indexed rows and columns by concentration position, which is only
!! right when the water's species are in stoichiometric order, not after a primary swap).
!!
!! Constraint types:
!!  icon 1 (concentration): c_k - ctot_k
!!  icon 2 (aqueous component): sum_i U_aq(k,i) c_i - ctot_k, over the aqueous species (primaries, then
!!         aqueous secondaries) with U_aq = [I | Se_1_star(aq)^T] (comp_mat_cst_act)
!!  icon 3 (activity): gamma_k c_k - ctot_k
!!  icon 4 (phase): sum_j Se(r,j) log10(a_j) + log10(ctot_k) - logK_r for the constraining reaction r
!! d(log10 gamma)/d(log10 c) enters through log_Jacobian_act_coeffs (rows by concentration index; columns
!! = [concentration index of the primaries | num_prim_species + aqueous secondary i], as built by
!! initialise_conc_anal) and dc2aq_dc1 (rows = aqueous secondaries, columns = primaries, stoichiometric).
subroutine compute_res_Jac_res_anal(this,indices_icon,n_icon,indices_constrains,ctot,dc2aq_dc1,log_Jacobian_act_coeffs,res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use arrays_m, only: int_array_c
    use array_ops_m, only: get_ind_int !> map a chemical-system eq. reaction index to its per-water logK position
    implicit none
    !> Pre-process
    class(aqueous_chemistry_c) :: this
    class(int_array_c), intent(in) :: indices_icon
    integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> (l,1) = concentration index of the constrained species, (l,2) = constraining eq. reaction index in the chemical system
    real(kind=8), intent(in) :: ctot(:)
    real(kind=8), intent(in) :: dc2aq_dc1(:,:) !> Jacobian aqueous secondary-aqueous primary
    real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
    real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson

    integer(kind=4) :: i,j,k,jk,ind_cstr,n_p,n_p_aq,n_2_aq,r
    integer(kind=4) :: jw !> row/column of water when its activity is prescribed (icon 3); 0 otherwise
    integer(kind=4), allocatable :: p(:),s(:),row_of_conc(:)
    real(kind=8), allocatable :: U_aq(:,:),c_aq(:),dlg_dlogc1(:)
    real(kind=8) :: ln10,c_k,gam_k

    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    ln10=log(1d1)
    ind_cstr=0
    jw=0
    res=0d0
    Jac_res=0d0 !> Jacobian residue
    !> Species storage: primary j at p(j), aqueous secondary i at s(i); row_of_conc(k) = primary stored at k
    p=this%ind_prim_species(1:n_p_aq)
    s=this%ind_sec_species(1:n_2_aq)
    allocate(row_of_conc(size(this%concentrations)))
    row_of_conc=0
    do j=1,n_p_aq
        row_of_conc(p(j))=j
    end do
    !> Aqueous component matrix [I | Se_1_star(aq)^T] over [primaries | aqueous secondaries] (stoichiometric order)
    allocate(U_aq(n_p_aq,n_p_aq+n_2_aq),c_aq(n_p_aq+n_2_aq),dlg_dlogc1(n_p_aq))
    U_aq(:,1:n_p_aq)=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq,1:n_p_aq)
    U_aq(:,n_p_aq+1:n_p_aq+n_2_aq)=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:n_p_aq,n_p+1:n_p+n_2_aq)
    c_aq(1:n_p_aq)=this%concentrations(p)
    c_aq(n_p_aq+1:n_p_aq+n_2_aq)=this%concentrations(s)
    !> Concentrations: res_k = c_k - ctot_k
    do i=1,n_icon(1)
        k=indices_icon%cols(1)%col_1(i)
        jk=row_of_conc(k)
        res(jk)=this%concentrations(k)-ctot(k)
        Jac_res(jk,jk)=1d0
    end do
    !> Aqueous components: res_k = U_aq(k,:) c_aq - ctot_k; d/dc1_j = U_aq(k,j) + U_aq(k,sec) dc2aq_dc1(:,j)
    do i=1,n_icon(2)
        k=indices_icon%cols(2)%col_1(i)
        jk=row_of_conc(k)
        res(jk)=dot_product(U_aq(jk,:),c_aq)-ctot(k)
        do j=1,n_p_aq
            Jac_res(jk,j)=U_aq(jk,j)+dot_product(U_aq(jk,n_p_aq+1:n_p_aq+n_2_aq),dc2aq_dc1(:,j))
        end do
    end do
    !> Activities: res_k = gamma_k c_k - ctot_k;
    !> d/dc1_j = gamma_k c_k [d(log gamma_k)/d(log c1_j) + sum_i d(log gamma_k)/d(log c2_i) d(log c2_i)/d(log c1_j)] / c1_j + gamma_k delta_kj
    !> Water: the solvent. Its molality (c_w = 1000/18.015 mol/kg) is not an unknown of the speciation and its
    !> activity is given by the model, Raoult's law a_w = 1 - 0.018 sum_{i/=w} c_i (compute_act_water), a function
    !> of the dissolved solids only: "a_w = ctot" cannot be an equation for c_w (zero derivative with respect
    !> to c_w, singular Jacobian, and the constant offset a_w - ctot could never reach abs_tol). The icon 3 input
    !> for water (CHEPROO convention: the water activity is prescribed, not solved for) therefore gives the
    !> trivial row Delta c_w = 0, and the water column is removed from the other rows below (their mass action
    !> laws see log a_w as a constant, so the 1/(c_w ln10) chain-rule term of the phase rows does not apply).
    do i=1,n_icon(3)
        k=indices_icon%cols(3)%col_1(i)
        jk=row_of_conc(k)
        if (this%aq_phase%ind_wat>0 .and. k==this%indices_aq_species(this%aq_phase%ind_wat)) then
            jw=jk
            res(jk)=0d0
            Jac_res(jk,:)=0d0
            Jac_res(jk,jk)=1d0
            cycle
        end if
        gam_k=1d1**this%log_act_coeffs(k)
        c_k=this%concentrations(k)
        res(jk)=gam_k*c_k-ctot(k)
        call d_log_gamma_d_log_c1(k,dlg_dlogc1)
        do j=1,n_p_aq
            Jac_res(jk,j)=gam_k*c_k*dlg_dlogc1(j)/this%concentrations(p(j))
        end do
        Jac_res(jk,jk)=Jac_res(jk,jk)+gam_k
    end do
    !> Phases: res_k = Se(r,1:n_p) log10(a(primaries)) + log10(ctot_k) - logK_r;
    !> d/dc1_j = [Se(r,j) + sum_l Se(r,l) d(log gamma_l)/d(log c1_j)] / (c1_j ln10)
    do i=1,n_icon(4)
        ind_cstr=ind_cstr+1 !> index constrain
        k=indices_icon%cols(4)%col_1(i)
        jk=row_of_conc(k)
        r=indices_constrains(ind_cstr,2)
        associate(Se_r => this%solid_chemistry%reactive_zone%chem_syst%Se(r,1:n_p_aq))
        res(jk)=dot_product(Se_r,log10(this%activities(p)))+log10(ctot(k))-&
            this%logK(get_ind_int(r,this%solid_chemistry%reactive_zone%ind_eq_reacts))
        Jac_res(jk,:)=Se_r
        do j=1,n_p_aq
            if (Se_r(j)==0d0) cycle
            call d_log_gamma_d_log_c1(p(j),dlg_dlogc1)
            Jac_res(jk,:)=Jac_res(jk,:)+Se_r(j)*dlg_dlogc1
        end do
        do j=1,n_p_aq
            Jac_res(jk,j)=Jac_res(jk,j)/(this%concentrations(p(j))*ln10)
        end do
        end associate
    end do
    !> Prescribed water activity: c_w is fixed, no row may depend on it (see the icon 3 block)
    if (jw>0) then
        Jac_res(:,jw)=0d0
        Jac_res(jw,jw)=1d0
    end if
    deallocate(p,s,row_of_conc,U_aq,c_aq,dlg_dlogc1)
contains
    !> Total derivative d(log10 gamma_k)/d(log10 c1_j) of the species stored at concentration index k with
    !> respect to every primary j (stoichiometric order): the direct term plus the one through the aqueous
    !> secondaries, d(log c2_i)/d(log c1_j) = dc2aq_dc1(i,j) c1_j / c2_i
    subroutine d_log_gamma_d_log_c1(k,dlg)
        integer(kind=4), intent(in) :: k
        real(kind=8), intent(out) :: dlg(:)
        integer(kind=4) :: jj,ii
        do jj=1,n_p_aq
            dlg(jj)=log_Jacobian_act_coeffs(k,p(jj))
            do ii=1,n_2_aq
                if (this%concentrations(s(ii))>0d0) dlg(jj)=dlg(jj)+log_Jacobian_act_coeffs(k,n_p+ii)*dc2aq_dc1(ii,jj)*&
                    this%concentrations(p(jj))/this%concentrations(s(ii))
            end do
        end do
    end subroutine d_log_gamma_d_log_c1
end subroutine
