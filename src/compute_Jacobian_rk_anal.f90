!> \file compute_Jacobian_rk_anal.f90
!> \brief Computes the Jacobian of kinetic reaction rates with respect to aqueous species concentrations using analytical derivatives
!> \details This subroutine calculates the Jacobian matrix ∂rₖ/∂c for all kinetic reactions in the system,
!> where rₖ are the kinetic reaction rates and c are the aqueous species concentrations. The Jacobian is
!> computed analytically (not via finite differences) for three types of kinetic reactions:
!> - Linear kinetic reactions (first-order or pseudo-first-order)
!> - Mineral kinetic reactions (dissolution/precipitation with TST-based rate laws)
!> - Redox kinetic reactions (Monod kinetics for microbial processes)
!>
!> The Jacobian is essential for implicit or semi-implicit time integration methods and Newton-Raphson
!> solvers in reactive transport simulations.
!>
!> Mathematical formulation:
!> \f[
!> \mathbf{J}_{rk} = \frac{\partial \mathbf{r}_k}{\partial \mathbf{c}} = 
!> \begin{bmatrix}
!> \frac{\partial r_{k,1}}{\partial c_1} & \cdots & \frac{\partial r_{k,1}}{\partial c_n} \\
!> \vdots & \ddots & \vdots \\
!> \frac{\partial r_{k,m}}{\partial c_1} & \cdots & \frac{\partial r_{k,m}}{\partial c_n}
!> \end{bmatrix}
!> \f]
!> where m is the number of kinetic reactions and n is the number of aqueous species.
!>
!> \param[in] this Aqueous chemistry object containing solid chemistry, reactive zone, chemical system, and kinetic reaction data
!> \param[out] drk_dc Jacobian matrix: ∂rₖ/∂c (must be already allocated with dimensions num_kin_reactions × num_aq_species) [T⁻¹]

subroutine compute_Jacobian_rk_anal(this,drk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !> Import aqueous chemistry class containing all kinetic reaction information
    implicit none !> Enforce explicit variable declarations
    class(aqueous_chemistry_c), intent(in) :: this !> Aqueous chemistry object containing solid chemistry, reactive zone, chemical system, and kinetic reactions [-]
    real(kind=8), intent(out) :: drk_dc(:,:) !> Jacobian matrix: ∂rₖ/∂c (must be already allocated, dimensions: num_kin_reactions × num_aq_species) [T⁻¹]
    
    integer(kind=4) :: i !> Loop counter for kinetic reactions of each type [-]
    integer(kind=4) :: k !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: n !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: m !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: l !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: rk_ind !> Kinetic reaction index (not used in current implementation) [-]
    integer(kind=4) :: cat_ind !> Catalyser index (not used in current implementation) [-]
    integer(kind=4) :: p !> Loop counter (not used in current implementation) [-]
    integer(kind=4) :: inh_ind !> Inhibitor index (not used in current implementation) [-]
    integer(kind=4) :: DOC_ind !> Dissolved organic carbon index (not used in current implementation) [-]
    integer(kind=4) :: num_kin !> Running count of total kinetic reactions processed (linear + mineral + redox) [-]
    integer(kind=4), allocatable :: indices(:) !> Indices of aqueous species participating in reaction (not used in current implementation) [-]
    real(kind=8), allocatable :: act_cat(:) !> Activities of catalyser species (not used in current implementation) [-]
    real(kind=8), allocatable :: drk_dc_loc(:) !> Local gradient vector for current kinetic reaction: ∂rₖ,ᵢ/∂c [T⁻¹]
    real(kind=8) :: saturation !> Saturation index for mineral kinetic reactions: Ω = Q/K_eq [-]

    drk_dc=0d0 !> Initialize Jacobian matrix to zero (default for species not participating in reactions)
    num_kin=0 !> Initialize counter for total number of kinetic reactions processed
!> Linear kinetic reactions: first-order or pseudo-first-order reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !> Loop over all linear kinetic reactions in chemical system
        num_kin=num_kin+1 !> Increment total kinetic reaction counter
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc(num_kin,:)) !> Compute gradient for i-th linear kinetic reaction: ∂rₖ/∂c (typically sparse, only non-zero for reactants/products)
    end do !> End loop over linear kinetic reactions
    !num_kin=num_kin+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !> Commented: Alternative way to increment counter after loop
!> Mineral kinetic reactions: dissolution/precipitation reactions with TST-based rate laws
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin !> Loop over all kinetic minerals in mineral zone
        num_kin=num_kin+1 !> Increment total kinetic reaction counter
        !indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%indices_aq_phase !> Commented: Extract indices of aqueous species in reaction
        allocate(drk_dc_loc(this%indices_rk%cols(num_kin)%dim)) !> Allocate local gradient vector with dimension = number of species participating in this reaction
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i)) !> Compute saturation index Ω = Q/K_eq for i-th kinetic mineral
        call this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Call mineral kinetic reaction gradient computation
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_drk_dc_mineral(this%concentrations(& !> Pass concentrations of aqueous species participating in reaction
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Extract species indices for this mineral reaction
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !> Map to aqueous phase indices
            this%activities(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Pass activities of aqueous species participating in reaction
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !> Map to aqueous phase indices for activities
            this%log_act_coeffs(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Pass log₁₀(activity coefficients) of aqueous species
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !> Map to aqueous phase indices for activity coefficients
            this%activities(this%indices_aq_species(& !> Pass activities of catalyser species (H⁺, OH⁻, etc.)
            this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%params%cat_indices)),saturation,& !> Include saturation index for TST rate law
            this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,drk_dc_loc) !> Pass reactive surface area and temperature, receive gradient in drk_dc_loc
        drk_dc(num_kin,this%indices_rk%cols(num_kin)%col_1)=drk_dc_loc !> Store local gradient in full Jacobian at appropriate columns (temporary workaround for sparse storage)
        deallocate(drk_dc_loc) !> Deallocate local gradient vector before next iteration
    end do !> End loop over mineral kinetic reactions
    !num_kin=num_kin+this%solid_chemistry%mineral_zone%num_minerals_kin !> Commented: Alternative way to increment counter after loop
!> Redox kinetic reactions: Monod kinetics for microbial processes (electron acceptor, electron donor, inhibition)
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts !> Loop over all redox kinetic reactions in chemical system
        !indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase !> Commented: Extract indices of aqueous species in reaction
        allocate(drk_dc_loc(this%indices_rk%cols(num_kin)%dim)) !> Allocate local gradient vector for this redox reaction
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(& !> Call Monod kinetics gradient computation
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),& !> Pass concentrations of species in Monod rate law (inhibitors, acceptor, donor)
            this%rk_new(num_kin-this%solid_chemistry%mineral_zone%num_minerals_kin),drk_dc_loc) !> Pass current reaction rate (offset by number of mineral reactions) and receive gradient
        drk_dc(num_kin,this%indices_rk%cols(num_kin)%col_1)=drk_dc_loc !> Store local gradient in full Jacobian at appropriate columns (temporary workaround for sparse storage)
        deallocate(drk_dc_loc) !> Deallocate local gradient vector before next iteration
    end do !> End loop over redox kinetic reactions
end subroutine