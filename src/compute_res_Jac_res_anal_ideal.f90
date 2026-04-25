!> \file compute_res_Jac_res_anal_ideal.f90
!! \brief Compute residual and Jacobian analytically for chemical equilibrium with ideal activity coefficients
!! \details This subroutine computes the residual vector and its Jacobian matrix for the Newton-Raphson
!! iteration used in chemical equilibrium calculations during initialization. It assumes ideal solution
!! behavior (activity coefficients γ = 1, so activities equal concentrations).
!!
!! The residual equations are constructed based on different constraint types (icon):
!! - **icon = 1**: Concentration constraints: r = c - c_total
!! - **icon = 2**: Aqueous component constraints: r = U·c - u_total
!! - **icon = 3**: Activity constraints: r = a - a_total
!! - **icon = 4**: Phase equilibrium constraints: r = log(K/a_phase) - Σ(S·log(a_primary))
!!
!! The Jacobian matrix ∂r/∂c is computed analytically:
!! \f[
!!   \frac{\partial r_i}{\partial c_j} = \begin{cases}
!!     1 & \text{for concentration constraints} \\
!!     U_{ij} + \sum_k U_{ik} \frac{\partial c_{2,k}}{\partial c_{1,j}} & \text{for component constraints} \\
!!     \gamma_i & \text{for activity constraints (γ=1 in ideal case)} \\
!!     \frac{S_{ij}}{c_j \ln(10)} & \text{for phase equilibrium constraints}
!!   \end{cases}
!! \f]
!!
!! **Assumptions:**
!! - Ideal solution behavior: γᵢ = 1 for all species
!! - Activities aᵢ = cᵢ (no ionic strength correction)
!! - Used only during initialization to obtain initial guess
!!
!! \param[in] this Aqueous chemistry object containing species, reactions, and current state
!! \param[in] indices_icon Indices of constraints for each icon type (1=conc, 2=comp, 3=act, 4=phase)
!! \param[in] n_icon Number of constraints for each icon type [4 integers]
!! \param[in] indices_constrains Indices linking constraints to chemical system reactions [n_cstr x 2]
!! \param[in] ctot Total concentration/activity values for each constraint [C or -]
!! \param[in] dc2aq_dc1 Jacobian of aqueous secondary species w.r.t. primary species ∂c₂/∂c₁ [-]
!! \param[out] res Residual vector for Newton-Raphson iteration [C or -]
!! \param[out] Jac_res Jacobian matrix of residual w.r.t. primary variables ∂r/∂c [-]
!!
!! \use aqueous_chemistry_m Provides aqueous_chemistry_c type for aqueous chemical system
!! \use arrays_m Provides int_array_c type for integer array storage

!> Computes residual and its Jacobian analytically in initialisation assuming ideal activities (γ = 1)
subroutine compute_res_Jac_res_anal_ideal(this,indices_icon,n_icon,indices_constrains,ctot,dc2aq_dc1,res,Jac_res)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    use arrays_m, only: int_array_c
    implicit none !> Enforce explicit variable declaration (no implicit typing)
    !> Input/output variables
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object containing species concentrations, activities, and chemical system
    class(int_array_c), intent(in) :: indices_icon !> Indices of constraints for each icon type (cols(1)=conc, cols(2)=comp, cols(3)=act, cols(4)=phase)
    integer(kind=4), intent(in) :: n_icon(:) !> Number of constraints for each icon type: [n_conc, n_comp, n_act, n_phase] [-]
    integer(kind=4), intent(in) :: indices_constrains(:,:) !> Indices linking constraints to reactions in chemical system [n_constraints x 2]
    real(kind=8), intent(in) :: ctot(:) !< Total concentration or activity values for constraints (must be ordered as in input file) [C or -]
    real(kind=8), intent(in) :: dc2aq_dc1(:,:) !> Jacobian of aqueous secondary species w.r.t. primary species: ∂c₂ᵃᵍ/∂c₁ [-]
    real(kind=8), intent(out) :: res(:) !> Residual vector for Newton-Raphson iteration: r = f(c) - target [C or -]
    real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian matrix of residual w.r.t. primary variables: ∂r/∂c [-]
    
    real(kind=8), allocatable :: c1(:),c2(:),c2v_init(:),c1k(:),c2k(:),log_c2k(:),log_c2(:),out_prod(:,:) !> Concentration vectors (not used in this routine, declared for compatibility)
    real(kind=8), allocatable :: Delta_c1(:) !> Newton increment: Δc₁ = c₁^(i+1) - c₁^i (not used here) [C]
    real(kind=8), allocatable :: abs_tol_res(:) !> Absolute tolerances for residuals in Newton-Raphson (not used here) [C or -]
    real(kind=8), allocatable :: mat_lin_syst(:,:),U_aq(:,:),K(:),z2(:) !> Auxiliary matrices: U_aq = component matrix for aqueous species [-]
    integer(kind=4) :: i,j,ind_eqn,niter_Picard,ind_cstr !> Loop indices (i,j), equation index, Picard iterations (unused), constraint counter
    integer(kind=4), allocatable :: ind_aq_comp(:),cols(:),ind_aq_species(:),counters(:) !> Index arrays (not used in this routine)
    logical :: flag_gas,flag_min,flag_wat !> Flags for gas, mineral, water phases (not used in this routine)
    
    ind_cstr=0 !> Initialize constraint counter to zero (will be incremented for phase equilibrium constraints) [-]
    Jac_res=0d0 !> Initialize Jacobian matrix to zero (will be filled selectively for each constraint type) [-]
    allocate(U_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,& !> Allocate component matrix for aqueous species: rows = num primary species, cols = total aqueous species
        this%aq_phase%num_species)) !> Total number of aqueous species (primary + secondary)
    U_aq(:,1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)=& !> Copy primary species block: U_aq for primary = identity portion of full component matrix [-]
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:& !> Extract rows 1 to num_aq_prim_species
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,1:& !> Extract columns 1 to num_aq_prim_species (identity block)
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) !> This forms the identity submatrix for primary species
    U_aq(:,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)=& !> Copy secondary species block: U_aq for secondary = stoichiometric coefficients [-]
        this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act(1:& !> Extract rows 1 to num_aq_prim_species
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species,& !> For all primary species
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Start after all primary species (including non-aqueous)
        this%aq_phase%num_species) !> End at total number of aqueous species (secondary aqueous block)
    !> Icon type 1: Concentration constraints
    !print *, indices_icon%cols(1)%col_1 !> Debugging output: display indices of concentration constraints (commented out)
    do i=1,n_icon(1) !> Loop over all concentration constraints (icon = 1)
        res(indices_icon%cols(1)%col_1(i))=this%concentrations(indices_icon%cols(1)%col_1(i))-ctot(indices_icon%cols(1)%col_1(i)) !> Residual: r = c - c_total (difference between current and target concentration) [C]
        Jac_res(indices_icon%cols(1)%col_1(i),indices_icon%cols(1)%col_1(i))=1d0 !> Jacobian: ∂r/∂c = 1 (diagonal entry, direct dependence) [-]
    end do
    !> Icon type 2: Aqueous component constraints
    !print *, this%ind_prim_species !> Debugging output: display indices of primary species (commented out)
    !print *, this%ind_sec_species !> Debugging output: display indices of secondary species (commented out)
    do i=1,n_icon(2) !> Loop over all aqueous component constraints (icon = 2)
        res(indices_icon%cols(2)%col_1(i))=dot_product(U_aq(indices_icon%cols(2)%col_1(i),:),this%concentrations(& !> Residual: r = U·c - u_total (component balance equation) [C]
            [this%ind_prim_species,this%ind_sec_species]))-& !> Concatenate primary and secondary species concentrations (chapuza = workaround)
            ctot(indices_icon%cols(2)%col_1(i)) !> Subtract total component concentration target [C]
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> Loop over all aqueous primary species to compute Jacobian entries
            Jac_res(indices_icon%cols(2)%col_1(i),this%ind_prim_species(j))=U_aq(indices_icon%cols(2)%col_1(i),j)+dot_product(U_aq(& !> Jacobian: ∂r/∂c₁ = U₁ⱼ + Σₖ U₂ₖ·∂c₂ₖ/∂c₁ⱼ (chain rule for secondary species dependence) [-]
                indices_icon%cols(2)%col_1(i),this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:& !> Sum over secondary species contributions
                this%aq_phase%num_species),dc2aq_dc1(:,j)) !> Multiply component coefficients by Jacobian of secondary w.r.t. primary species
        end do
    end do
    !> Icon type 3: Activity constraints
    do i=1,n_icon(3) !> Loop over all activity constraints (icon = 3)
        res(indices_icon%cols(3)%col_1(i))=this%activities(indices_icon%cols(3)%col_1(i))-ctot(indices_icon%cols(3)%col_1(i)) !> Residual: r = a - a_total (difference between current and target activity) [-]
        Jac_res(indices_icon%cols(3)%col_1(i),indices_icon%cols(3)%col_1(i))=10**this%log_act_coeffs(indices_icon%cols(3)%col_1(i)) !> Jacobian: ∂a/∂c = γ = 10^(log₁₀γ) (activity coefficient, equals 1 in ideal case) [-]
    end do
    !> Icon type 4: Phase equilibrium constraints
    !print *, indices_icon%cols(4)%col_1 !> Debugging output: display indices of phase equilibrium constraints (commented out)
    do i=1,n_icon(4) !> Loop over all phase equilibrium constraints (icon = 4, typically mineral or gas phases)
        ind_cstr=ind_cstr+1 !> Increment constraint counter (maps icon constraint to chemical system reaction index) [-]
        !print *, ctot(indices_icon%cols(4)%col_1(i)) !> Debugging output: display total concentration/activity for this phase constraint
        !print *, this%solid_chemistry%reactive_zone%chem_syst%Se(& !> Debugging output: display stoichiometric matrix row for this equilibrium reaction
        !    indices_constrains(ind_cstr,2),1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) !> Row from secondary species matrix (equilibrium stoichiometry)
        !print *, this%activities(this%ind_prim_species) !> Debugging output: display current activities of primary species
        !print *, this%concentrations(this%ind_prim_species) !> Debugging output: display current concentrations of primary species
        !print *, this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr,2))%eq_cst !> Debugging output: display equilibrium constant K for this reaction
        res(indices_icon%cols(4)%col_1(i))=dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(& !> Residual: r = Σ(Sᵢⱼ·log₁₀(aⱼ)) + log₁₀(a_phase) - log₁₀(K) (mass action law) [-]
            indices_constrains(ind_cstr,2),1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),& !> Stoichiometric coefficients for primary species in this reaction
            log10(this%activities(this%ind_prim_species)))+log10(ctot(indices_icon%cols(4)%col_1(i)))-& !> Add log₁₀ of phase activity (saturation index contribution)
            log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr,2))%eq_cst) !> Subtract log₁₀ of equilibrium constant K
        !print *, dot_product(this%solid_chemistry%reactive_zone%chem_syst%Se(& !> Debugging output: display first term of residual (stoichiometric contribution)
        !    indices_constrains(ind_cstr,2),1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species),& !> Primary species stoichiometry
        !    log10(this%activities(this%ind_prim_species))) !> Primary species log activities
        !print *, log10(ctot(indices_icon%cols(4)%col_1(i))) !> Debugging output: display second term of residual (log of phase activity)
        !print *, log10(this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr,2))%eq_cst) !> Debugging output: display third term of residual (log of equilibrium constant)
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> Loop over all aqueous primary species to compute Jacobian entries
            !print *, this%ind_prim_species(j) !> Debugging output: display index of current primary species
            Jac_res(indices_icon%cols(4)%col_1(i),this%ind_prim_species(j))=this%solid_chemistry%reactive_zone%chem_syst%Se(& !> Jacobian: ∂r/∂cⱼ = Sᵢⱼ/(cⱼ·ln(10)) (derivative of log term w.r.t. concentration) [1/C]
                indices_constrains(ind_cstr,2),j)/(this%concentrations(this%ind_prim_species(j))*log(1d1)) !> Divide stoichiometric coefficient by (concentration × ln(10)) for log derivative
        end do
        ! !> [DEBUG icon=4] dump residual / equilibrium info / Jacobian row for this constraint
        ! print '(A,I0,A,I0,A,A,A,ES13.5,A,ES13.5,A,ES13.5)', &
        !     '[DBG icon4] i=', i, &
        !     ' eq_react_idx=', indices_constrains(ind_cstr,2), &
        !     ' constrain_sp=', trim(this%aq_phase%aq_species(indices_icon%cols(4)%col_1(i))%name), &
        !     ' res=', res(indices_icon%cols(4)%col_1(i)), &
        !     ' ctot=', ctot(indices_icon%cols(4)%col_1(i)), &
        !     ' eq_cst=', this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(indices_constrains(ind_cstr,2))%eq_cst
        ! print '(A)', '   Se row (prim cols):'
        ! do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        !     print '(A,I0,A,A,A,ES13.5,A,ES13.5,A,ES13.5)', &
        !         '     j=', j, ' prim_sp=', trim(this%aq_phase%aq_species(this%ind_prim_species(j))%name), &
        !         ' Se=', this%solid_chemistry%reactive_zone%chem_syst%Se(indices_constrains(ind_cstr,2),j), &
        !         ' c=', this%concentrations(this%ind_prim_species(j)), &
        !         ' Jac=', Jac_res(indices_icon%cols(4)%col_1(i),this%ind_prim_species(j))
        ! end do
    end do
end subroutine !> End of compute_res_Jac_res_anal_ideal subroutine