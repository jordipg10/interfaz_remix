!> @file compute_rk_Jac_rk_anal.f90
!> @brief Computes kinetic reaction rates and their Jacobian matrix analytically
!> @details This subroutine computes reaction rates and their derivatives for all kinetic reactions
!> in an aqueous chemistry system, including linear kinetics, Monod kinetics, and mineral kinetics.
!> The analytical Jacobian drk_dc represents the sensitivity of each reaction rate to changes in 
!> aqueous concentrations, essential for implicit time-stepping methods (Newton-Raphson solvers).
!>
!> @par Algorithm Overview:
!> 1. Linear kinetic reactions: First-order reactions with constant rate coefficients
!> 2. Redox kinetic reactions: Monod-type rate laws with multiple reactant dependencies
!> 3. Mineral kinetic reactions: Surface-area normalized rate laws with saturation state dependencies
!>
!> @par Mathematical Context:
!> For each reaction k, computes:
!> - rk_new(k): Reaction rate [mol/(L·s)]
!> - drk_dc(k,j): Jacobian ∂rk/∂c_j [1/s]
!>
!> @see compute_rk_Jac_rk_incr_coeff For numerical Jacobian approximation
!> @see aqueous_chemistry_c For aqueous chemistry data structure
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

!> @brief Computes kinetic reaction rates and analytical Jacobian for aqueous chemistry
!> @param[in,out] this Aqueous chemistry object containing species concentrations, activities, and chemistry system
!> @param[out] rk_new Array of newly computed kinetic reaction rates [mol/(L·s)] (must be pre-allocated)
!> @param[out] drk_dc Jacobian matrix ∂rk_i/∂c_j [1/s] (must be pre-allocated, sparse format via indices)
subroutine compute_rk_Jac_rk_anal(this,rk_new,drk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !< Import aqueous chemistry class definition
    implicit none !< Enforce explicit variable declarations for type safety
    
!> @section args Subroutine Arguments
    class(aqueous_chemistry_c) :: this !< Aqueous chemistry object with species data and reaction system
    real(kind=8), intent(out) :: rk_new(:) !< Output array for kinetic reaction rates [mol/(L·s)] (must be allocated)
    real(kind=8), intent(out) :: drk_dc(:,:) !< Output Jacobian matrix ∂rk/∂c [1/s] (must be allocated, sparse storage)
    
!> @section vars Local Variables
    integer(kind=4) :: i !< Loop counter for iterating over reactions within each reaction type
    integer(kind=4) :: n !< Unused integer variable (legacy code, can be removed)
    integer(kind=4) :: niter !< Unused integer variable (legacy code, can be removed)
    integer(kind=4) :: rk_ind !< Unused integer variable (legacy code, can be removed)
    integer(kind=4) :: l !< Unused integer variable (legacy code, can be removed)
    integer(kind=4) :: index !< Unused integer variable (legacy code, can be removed)
    integer(kind=4) :: num_rk !< Running counter for total number of kinetic reactions processed so far
    integer(kind=4), allocatable :: indices(:) !< Unused allocatable array (legacy code, can be removed)
    real(kind=8), allocatable :: drk_dc_loc(:) !< Local temporary array for storing Jacobian derivatives of single reaction
    real(kind=8) :: saturation !< Saturation state Ω = IAP/K_eq for mineral dissolution/precipitation kinetics

    !> Debug: print entry information
    ! print *, "DEBUG compute_rk_Jac_rk_anal: entering"
    ! print *, "DEBUG compute_rk_Jac_rk_anal: size(rk_new) =", size(rk_new)
    ! print *, "DEBUG compute_rk_Jac_rk_anal: shape(drk_dc) =", shape(drk_dc)
    ! print *, "DEBUG compute_rk_Jac_rk_anal: num_lin_kin_reacts =", &
    !     this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
    ! print *, "DEBUG compute_rk_Jac_rk_anal: num_redox_kin_reacts =", &
    !     this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    ! print *, "DEBUG compute_rk_Jac_rk_anal: num_minerals_kin =", &
    !     this%solid_chemistry%mineral_zone%num_minerals_kin

    num_rk=0 !< Initialize counter for kinetic reactions (starts at zero before processing any reactions)
    
    !> @note Legacy code: Old rate updates are commented out (may be handled elsewhere in calling code)
    !!> Update old kinetic reaction rates
    !call this%update_rk_old() !> we update old values of kinetic reaction rates
    !call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
    
!> @subsection lin_kin Linear Kinetic Reactions
!> @details Process all linear kinetic reactions (first-order rate laws: rk = k * c)
!> Linear kinetics: rk = k * c_i, where k is rate constant and c_i is reactant concentration
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- BEGIN linear kinetic reactions ---"
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !< Loop over all linear kinetic reactions
        num_rk=num_rk+1 !< Increment global reaction counter (tracks position in rk_new and drk_dc arrays)
        ! print *, "DEBUG compute_rk_Jac_rk_anal: linear reaction i =", i, " num_rk =", num_rk
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   indices_aq_phase(1) =", &
        !     this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   concentration =", &
        !     this%concentrations(this%indices_aq_species(&
        !     this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)))
        !index=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)
        allocate(drk_dc_loc(1)) !< Allocate temporary array for single derivative (linear kinetics has one reactant)
        !> Compute linear reaction rate: rk = k * c (first-order kinetics)
        !> Extract concentration of the single reactant species from full concentration array
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_rk_lin(&
            this%concentrations(this%indices_aq_species(&
            this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1))),& !< Reactant concentration [mol/L]
            this%rk_new(i)) !< Store rate in chemistry object
        rk_new(i)=this%rk_new(i) !< Copy rate to output array (temporary workaround - "chapuza")
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   rk_new(", i, ") =", rk_new(i)
        allocate(drk_dc_loc(size(this%indices_rk%cols(num_rk)%col_1))) !< Reallocate for sparse Jacobian column
        !> Compute Jacobian derivative drk/dc = k (constant for linear kinetics)
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc_loc)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   drk_dc_loc =", drk_dc_loc
        !> Store Jacobian entry in sparse matrix format using column indices
        drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1) = drk_dc_loc !< Assign derivatives to sparse column
        deallocate(drk_dc_loc) !< Free temporary memory for next iteration
        ! drk_dc(i,this%indices_aq_species(index))=drk_dc_loc(1) !> chapuza
        !deallocate(drk_dc_loc)
    end do !< End of linear kinetic reactions loop
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- END linear kinetic reactions, num_rk =", num_rk
    !num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !< Legacy: redundant increment
    
!> @subsection redox_kin Redox Kinetic Reactions (Monod-Type Rate Laws)
!> @details Process all redox kinetic reactions using Monod/Michaelis-Menten kinetics
!> Monod kinetics: rk = r_max * [Π(c_i/(K_i + c_i))] - accounts for multiple limiting substrates
!> @par Mathematical Form:
!> rk = V_max * c_donor/(K_donor + c_donor) * c_acceptor/(K_acceptor + c_acceptor) * biomass
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- BEGIN redox kinetic reactions ---"
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts !< Loop over all redox kinetic reactions
        num_rk=num_rk+1 !< Increment global reaction counter
        ! print *, "DEBUG compute_rk_Jac_rk_anal: redox reaction i =", i, " num_rk =", num_rk
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   indices_aq_phase =", &
        !     this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   concentrations =", &
        !     this%concentrations(this%indices_aq_species(&
        !     this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase))
        ! indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        allocate(drk_dc_loc(this%indices_rk%cols(num_rk)%dim)) !< Allocate for multiple derivatives (Monod depends on multiple species)
        !> @note Legacy code: Old implementation split rate and Jacobian computation (now combined for efficiency)
        !call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(this%concentrations(&
        !    this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),&
        !    this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin)) !> chapuza
        !rk(num_rk)=this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin) !> chapuza
        !call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(&
        !    this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),&
        !    this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin),drk_dc_loc)
        !> Compute Monod reaction rate and Jacobian simultaneously (optimized approach)
        !> Extract concentrations of all species involved in the redox reaction (electron donors/acceptors, biomass)
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_drk_dc_Monod(&
            this%concentrations(this%indices_aq_species(&
            this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),& !< Substrate/biomass concentrations [mol/L]
            rk_new(num_rk),& !< Output: reaction rate [mol/(L·s)]
            drk_dc_loc) !< Output: Jacobian derivatives [1/s]
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   rk_new(", num_rk, ") =", rk_new(num_rk)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   drk_dc_loc =", drk_dc_loc
        !> Apply minimum rate threshold to prevent numerical issues in stiff solvers
        if (abs(rk_new(num_rk))<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then !< Check if rate is below convergence tolerance squared
            ! print *, "DEBUG compute_rk_Jac_rk_anal:   rate below threshold, clamping to abs_tol^2 =", &
            !     this%solid_chemistry%reactive_zone%CV_params%abs_tol**2
            this%rk_new(num_rk)=&
                this%solid_chemistry%reactive_zone%CV_params%abs_tol**2 !< Set minimum rate floor (temporary workaround - "chapuza")
        else
            this%rk_new(num_rk)=rk_new(num_rk) !< Use computed rate (temporary workaround - "chapuza")
        end if
        drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1)=drk_dc_loc !< Store Jacobian in sparse format (temporary workaround - "chapuza")
        deallocate(drk_dc_loc) !< Free temporary memory for next iteration
    end do !< End of redox kinetic reactions loop
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- END redox kinetic reactions, num_rk =", num_rk
    
!> @subsection min_kin Mineral Kinetic Reactions (Dissolution/Precipitation)
!> @details Process all mineral kinetic reactions with surface-area normalized rate laws
!> Mineral kinetics typically follow transition state theory (TST):
!> rk = A * k * Π(a_i^n_i) * f(Ω) where:
!>   - A = reactive surface area [m²/L]
!>   - k = rate constant [mol/(m²·s)]
!>   - a_i = activity of catalyzing/inhibiting species
!>   - n_i = reaction order
!>   - f(Ω) = saturation state function, e.g., (1 - Ω) or (1 - Ω^n)
!> @par Saturation State:
!> Ω = IAP/K_eq where IAP = ion activity product, K_eq = equilibrium constant
!> Ω < 1: undersaturated (dissolution), Ω > 1: supersaturated (precipitation)
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- BEGIN mineral kinetic reactions ---"
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin !< Loop over all kinetic minerals in the system
        ! print *, "DEBUG compute_rk_Jac_rk_anal: mineral reaction i =", i, " num_rk+i =", num_rk+i
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   ind_min_chem_syst(i) =", &
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   indices_aq_phase =", &
        !     this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   cat_indices =", &
        !     this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   react_surface(i) =", &
        !     this%solid_chemistry%react_surfaces(i)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   temp =", this%solid_chemistry%temp
        !num_rk=num_rk+1 !< Legacy: counter increment now done after loop
        !indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
        !    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase
        allocate(drk_dc_loc(this%indices_rk%cols(num_rk+i)%dim)) !< Allocate for Jacobian column (mineral kinetics depends on multiple ions)
        !> Compute saturation state Ω = IAP/K_eq for this mineral
        !> Determines whether mineral dissolves (Ω<1) or precipitates (Ω>1)
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   saturation (Omega) =", saturation
        !> Compute mineral reaction rate using TST-based rate law
        !> Rate depends on: catalyzing ion activities, saturation state, reactive surface area, temperature
        call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_rk_mineral(&
            this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),& !< Activities of catalyzing ions (e.g., H+ for acid dissolution)
            saturation,& !< Saturation state Ω = IAP/K_eq [-]
            this%solid_chemistry%react_surfaces(i),& !< Reactive surface area [m²/L]
            this%solid_chemistry%temp,& !< Temperature [K] for Arrhenius temperature dependence
            this%solid_chemistry%rk_new(i)) !< Output: mineral reaction rate [mol/(L·s)]
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   rk_mineral (before time scaling) =", this%solid_chemistry%rk_new(i)
        !> chapuza Paulina
        ! print *, this%activities(this%indices_aq_species(&
        !     this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase(&
        !     this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk)))
        this%solid_chemistry%rk_new(i)=this%solid_chemistry%rk_new(i)*this%activities(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk)))*&
            this%solid_chemistry%time**(&
            -this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%m) !< Update rate in chemistry object (temporary workaround - "chapuza")
        rk_new(num_rk+i)=this%solid_chemistry%rk_new(i) !< Copy rate to output array (temporary workaround - "chapuza")
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   rk_new(", num_rk+i, ") (after time scaling) =", rk_new(num_rk+i)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   time =", this%solid_chemistry%time, &
        !     " m =", this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !     this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%m
        !> Compute Jacobian matrix drk/dc for mineral reaction
        !> Derivatives account for:
        !>   1. Effect of concentration on activity (via activity coefficients)
        !>   2. Effect of activity on rate (via catalysis terms and saturation function)
        !> Chain rule: drk/dc = drk/da * da/dc = drk/da * (γ + c * dγ/dc)
        call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_drk_dc_mineral(&
            this%concentrations(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !< Concentrations of all species involved [mol/L]
            this%activities(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !< Activities of all species involved [-]
            this%log_act_coeffs(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !< Log activity coefficients log(γ) for chain rule derivatives
            this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),& !< Activities of catalyzing ions
            saturation,& !< Saturation state Ω (precomputed above)
            this%solid_chemistry%react_surfaces(i),& !< Reactive surface area [m²/L]
            this%solid_chemistry%temp,& !< Temperature [K]
            drk_dc_loc) !< Output: Jacobian derivatives drk/dc [1/s]
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   drk_dc_loc (mineral) (before time scaling) =", drk_dc_loc
        drk_dc_loc(1:this%indices_rk%cols(num_rk+i)%dim-1)=drk_dc_loc(1:this%indices_rk%cols(num_rk+i)%dim-1)*&
            this%activities(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk)))*&
            this%solid_chemistry%time**(&
            -this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%m)
        drk_dc_loc(this%indices_rk%cols(num_rk+i)%dim)=-this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%k(1)*&
            this%solid_chemistry%time**(&
            -this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%m)
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   drk_dc_loc (mineral) (after time scaling) =", drk_dc_loc
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   indices_rk (mineral) =", &
        !     this%indices_rk%cols(num_rk+i)%col_1
        drk_dc(num_rk+i,this%indices_rk%cols(num_rk+i)%col_1)=drk_dc_loc !< Store Jacobian in sparse format (temporary workaround - "chapuza")
        ! print *, "DEBUG compute_rk_Jac_rk_anal:   drk_dc (mineral) =", drk_dc
        deallocate(drk_dc_loc) !< Free temporary memory for next iteration
    end do !< End of mineral kinetic reactions loop
    num_rk=num_rk+this%solid_chemistry%mineral_zone%num_minerals_kin !< Update total reaction counter with number of mineral reactions
    ! print *, "DEBUG compute_rk_Jac_rk_anal: --- END mineral kinetic reactions, num_rk =", num_rk
    ! print *, "DEBUG compute_rk_Jac_rk_anal: completed successfully, total rk_new =", rk_new(1:num_rk)
end subroutine !< End of compute_rk_Jac_rk_anal subroutine