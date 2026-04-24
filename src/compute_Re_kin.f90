!> Computes equilibrium reaction amounts from secondary variable activity species concentrations using linear least squares
!! Reaction rates are expressed per unit volume of water
subroutine compute_Re_kin(this,c2v_hat,Delta_t,lambda_r)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at time step k+1
    real(kind=8), intent(in) :: c2v_hat(:) !> concentrations secondary variable activity species after mixing at time step k
    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
    !real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
    real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio associated to "this" argument [-]
!> Variables
    real(kind=8), allocatable :: A(:,:),b(:),c2v(:),Re(:),Sk_nc(:,:),Rk(:),Se2v(:,:),Sk2v(:,:)
    integer(kind=4) :: err
!> Pre-process
    allocate(Re(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    allocate(A(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions,&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
    allocate(b(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
!> Process
    !> Linear least squares
    c2v=this%get_c2v()
    Rk=this%get_Rk() !> we get kinetic reaction amounts at time step k+1
    Sk_nc=this%get_Sk_nc()
    Sk2v=Sk_nc(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:)
    Se2v=this%solid_chemistry%reactive_zone%get_Se2v_react_zone() !> Retrieve stoichiometric matrix for kinetic reactions (columns for non-constant activity species) - retrieved but not used in current implementation [-]
    A=lambda_r*matmul(Se2v,transpose(Se2v))
    b=matmul(Se2v,c2v-c2v_hat-lambda_r*matmul(transpose(Sk2v),Rk)) !> b=Delta_t*(Re+matmul(transpose(Sk_nc),theta*rk+(1d0-theta)*rk_old))
    if (inf_norm_vec_real(b)<this%solid_chemistry%reactive_zone%CV_params%zero) then
        Re=0d0
    else
        call LU_lin_syst(A,b,this%solid_chemistry%reactive_zone%CV_params%zero,Re) !> linear system solver
    end if
    this%solid_chemistry%Re(1:this%solid_chemistry%reactive_zone%num_minerals_cst_act)=Re(1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act)!/Delta_t !> Re=Re/Delta_t
    this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals_cst_act+1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act+this%solid_chemistry%reactive_zone%num_minerals_var_act)=&
        Re(this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species-this%aq_phase%wat_flag+&
        this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species-&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)!/Delta_t !> Re=Re/Delta_t
    this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals+1:this%solid_chemistry%reactive_zone%num_minerals+&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=&
        Re(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species-&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species)!/Delta_t !> Re=Re/Delta_t
    call this%solid_chemistry%set_re_mean(Delta_t) !> Update mean equilibrium reaction rates in solid chemistry object
    if (associated(this%gas_chemistry)) then
        this%gas_chemistry%Re(1:this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)=Re(&
            this%solid_chemistry%reactive_zone%num_minerals_cst_act+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species-this%aq_phase%wat_flag) !> Re_j=Re/Delta_t
        this%gas_chemistry%Re(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=Re(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)!/Delta_t !> Re=Re/Delta_t
        call this%gas_chemistry%set_re_mean(Delta_t) !> Update mean equilibrium reaction rates in gas chemistry object
    end if
    this%Re=Re(this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+&
        this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts)!/Delta_t !> Re_j=Re/Delta_t
    call this%set_re_mean(Delta_t) !> Update mean equilibrium reaction rates in aqueous chemistry object
    deallocate(A,b,c2v,Re,Sk_nc,Se2v)
end subroutine