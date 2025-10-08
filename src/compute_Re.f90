!> Computes equilibrium reaction amounts from secondary variable activity species concentrations using linear least squares
!! Reaction rates are expressed per unit volume of water
!! We DO NOT apply lumping to the mixing ratios of reaction rates
subroutine compute_Re(this,c2nc_tilde,Delta_t,theta)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at time step k+1
    real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
    real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
    !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions to secondary variable activity species after mixing
!> Variables
    real(kind=8), allocatable :: A(:,:),b(:),c2nc(:),Re(:),Sk_nc(:,:),rk_new(:),rk_old(:)
    integer(kind=4) :: err
!> Pre-process
    allocate(Re(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)) !> re_mean=Delta_t*re_mean
!> Process
    !> Linear least squares
    c2nc=this%get_c2nc()
    ! rk=this%get_rk() !> we get kinetic reaction rates (chapuza)
    ! rk_old=this%get_rk_old() !> we get old kinetic reaction rates (chapuza)
    Sk_nc=this%get_Sk_nc()
    A=matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),transpose(&
        this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)))
    b=matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),c2nc-c2nc_tilde)
        ! matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,&
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),Delta_t*(rk_tilde+matmul(transpose(&
        ! Sk_nc(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)),&
        ! theta*rk+(1d0-theta)*rk_old))) !> b=Delta_t*(re_mean+matmul(transpose(Sk_nc),theta*rk+(1d0-theta)*rk_old))
    if (inf_norm_vec_real(b)<this%solid_chemistry%reactive_zone%CV_params%zero) then
        Re=0d0
    else
        call LU_lin_syst(A,b,this%solid_chemistry%reactive_zone%CV_params%zero,Re) !> linear system solver
    end if
    this%solid_chemistry%Re(1:this%solid_chemistry%reactive_zone%num_minerals_cst_act)=Re(1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act)
    this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals_cst_act+1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act+this%solid_chemistry%reactive_zone%num_minerals_var_act)=&
        Re(this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species-this%aq_phase%wat_flag+&
        this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species-&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
    this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals+1:this%solid_chemistry%reactive_zone%num_minerals+&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=&
        Re(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species-&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1:&
        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species)!/Delta_t !> re_mean=re_mean/Delta_t
    if (associated(this%gas_chemistry)) then
        this%gas_chemistry%Re(1:this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)=Re(&
            this%solid_chemistry%reactive_zone%num_minerals_cst_act+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species-this%aq_phase%wat_flag) !> re_mean_j=re_mean/Delta_t
        this%gas_chemistry%Re(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=Re(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)!/Delta_t !> re_mean=re_mean/Delta_t
    end if
    this%Re=Re(this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
        this%solid_chemistry%reactive_zone%num_minerals_cst_act+&
        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+&
        this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts)!/Delta_t !> re_mean_j=re_mean/Delta_t
end subroutine