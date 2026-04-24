!> \file compute_rk_new.f90
!> \brief Computes all kinetic reaction rates for an aqueous chemistry system
!> \details This subroutine calculates the kinetic reaction rates for all types of kinetic reactions
!> associated with an aqueous chemistry object. It sequentially computes rates for:
!> 1. Linear kinetic reactions (first-order or pseudo-first-order reactions)
!> 2. Redox kinetic reactions (microbial reactions using Monod kinetics)
!> 3. Mineral kinetic reactions (dissolution/precipitation using transition state theory)
!>
!> The subroutine updates the internal reaction rate arrays (this%rk_new, this%solid_chemistry%rk_new)
!> and returns the consolidated rates in the output array rk_new.
!>
!> The rates are stored in order:
!> - Indices 1 to num_lin_kin_reacts: Linear kinetic reactions
!> - Indices (num_lin_kin_reacts+1) to (num_lin_kin_reacts+num_redox_kin_reacts): Redox kinetic reactions
!> - Indices (num_lin_kin_reacts+num_redox_kin_reacts+1) to end: Mineral kinetic reactions
!>
!> \param[in,out] this Aqueous chemistry object containing species concentrations, activities, and reaction parameters [-]
!> \param[out] rk_new Array of all kinetic reaction rates (must be pre-allocated with size = total number of kinetic reactions) [M/(L³·T)]

subroutine compute_rk_new(this,rk_new)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !> Import aqueous chemistry class containing species, reactions, and solid chemistry
    implicit none !> Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object with concentrations, activities, solid chemistry, and reaction zone information [-]
    real(kind=8), intent(out) :: rk_new(:) !> Output array of kinetic reaction rates (must be pre-allocated) [M/(L³·T)]
!> Variables
    integer(kind=4) :: i !> Loop counter for iterating over reactions [-]
    integer(kind=4) :: n !> General loop counter (declared but not used in current implementation) [-]
    integer(kind=4) :: niter !> Iteration counter (declared but not used in current implementation) [-]
    integer(kind=4) :: rk_ind !> Reaction index (declared but not used in current implementation) [-]
    integer(kind=4) :: l !> Loop counter (declared but not used in current implementation) [-]
    integer(kind=4) :: index !> Species index (declared but not used in current implementation) [-]
    integer(kind=4) :: num_rk !> Running counter for total number of kinetic reactions processed so far [-]
    integer(kind=4), allocatable :: indices(:) !> Array of aqueous phase species indices for a reaction (declared but not used in current implementation) [-]
    real(kind=8), allocatable :: drk_dc_loc(:) !> Local array for storing reaction rate gradients (declared but not used in current implementation) [L³/M]
    real(kind=8) :: saturation !> Mineral saturation index (ratio of ion activity product to equilibrium constant) [-]

    num_rk=0 !> Initialize counter for total kinetic reactions processed (starts at zero before processing any reactions) [-]
    
!!> Update old kinetic reaction rates (commented out - old values not updated in current implementation)
!    call this%update_rk_old() !> Would update old values of aqueous kinetic reaction rates for time integration
!    call this%solid_chemistry%update_rk_old() !> Would update old values of solid chemistry kinetic reaction rates for time integration
!> We compute linear kinetic reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !> Loop over all linear kinetic reactions in the chemical system
        num_rk=num_rk+1 !> Increment total reaction counter (linear reactions are indexed first: 1 to num_lin_kin_reacts) [-]
        !index=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1) !> Commented: would get aqueous species index for this reaction
        !allocate(drk_dc_loc(1)) !> Commented: would allocate gradient array for single-species linear reaction
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_rk_lin(& !> Compute linear kinetic reaction rate for reaction i
            this%concentrations(this%indices_aq_species(& !> Pass concentration of the aqueous species involved in this linear reaction [M/L³]
            this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1))),this%rk_new(i)) !> Store rate in this%rk_new(i) [M/(L³·T)]
        rk_new(i)=this%rk_new(i) !> Copy linear reaction rate to output array (workaround: chapuza) [M/(L³·T)]
        !allocate(drk_dc_loc(size(this%indices_rk%cols(num_rk)%col_1))) !> Commented: would allocate gradient array based on Jacobian structure
        !call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc_loc) !> Commented: would compute reaction rate gradient
        !drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1) = drk_dc_loc !> Commented: would store gradient in Jacobian matrix
        !deallocate(drk_dc_loc) !> Commented: would free gradient array
        ! drk_dc(i,this%indices_aq_species(index))=drk_dc_loc(1) !> Commented: alternative gradient storage (workaround: chapuza)
        !deallocate(drk_dc_loc) !> Commented: would free gradient array
    end do !> End linear kinetic reactions loop
    !num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts !> Commented: alternative way to update counter after loop
!> We compute redox reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts !> Loop over all redox kinetic reactions in the chemical system (microbial Monod reactions)
        num_rk=num_rk+1 !> Increment total reaction counter (redox reactions indexed after linear and mineral reactions) [-]
        ! indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase !> Commented: would get aqueous species indices for this redox reaction
        !allocate(drk_dc_loc(this%indices_rk%cols(num_rk)%dim)) !> Commented: would allocate gradient array for redox reaction
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(this%concentrations(& !> Compute Monod kinetic rate for redox reaction i
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),& !> Pass concentrations of inhibitors, electron acceptor, and electron donor [M/L³]
            this%rk_new(num_rk)) !> Store rate in this%rk_new (index adjusted for mineral offset: workaround: chapuza) [M/(L³·T)]
        rk_new(num_rk)=this%rk_new(num_rk) !> Copy Monod reaction rate to output array at position num_rk (workaround: chapuza) [M/(L³·T)]
        !call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(& !> Commented: would compute Monod rate gradient
        !    this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),& !> Commented: would pass concentrations for gradient
        !    this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin),drk_dc_loc) !> Commented: would pass current rate and output gradient
        !drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1)=drk_dc_loc !> Commented: would store gradient in Jacobian matrix (workaround: chapuza)
        !deallocate(drk_dc_loc) !> Commented: would free gradient array
    end do !> End Monod redox kinetic reactions loop
!> We compute mineral kinetic reaction rates
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin !> Loop over all kinetic minerals in the mineral zone (dissolution/precipitation reactions)
        !num_rk=num_rk+1 !> Commented: would increment counter inside loop (counter updated after loop instead)
        !indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Commented: would get aqueous species indices for this mineral reaction
        !    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase !> Commented: mapping from mineral zone index to chemical system index
        !allocate(drk_dc_loc(this%indices_rk%cols(num_rk+i)%dim)) !> Commented: would allocate gradient array for mineral reaction
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i)) !> Compute saturation index Ω = IAP/K_eq for kinetic mineral i (Ω>1: supersaturated, Ω<1: undersaturated) [-]
        call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(& !> Call mineral kinetic rate computation for mineral i
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_rk_mineral(& !> Access mineral reaction via chemical system index
            this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(& !> Pass activities of catalyser species (species affecting rate) [M/L³]
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation,& !> Pass saturation index Ω [-]
            this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,this%solid_chemistry%rk_new(i)) !> Pass reactive surface area [L²] and temperature [K], store rate in solid_chemistry%rk_new(i) [M/(L³·T)]
        rk_new(num_rk+i)=this%solid_chemistry%rk_new(i) !> Copy mineral reaction rate to output array at position num_rk+i (workaround: chapuza) [M/(L³·T)]
        !call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(& !> Commented: would compute mineral rate gradient
        !    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_drk_dc_mineral(& !> Commented: would call gradient computation
        !        this%concentrations(this%indices_aq_species(& !> Commented: would pass concentrations for gradient
        !        this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(& !> Commented: accessing reaction parameters
        !        this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),& !> Commented: species indices
        !        this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(& !> Commented: would pass activities
        !        this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation,& !> Commented: saturation for gradient
        !        this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,drk_dc_loc) !> Commented: surface, temperature, output gradient
        !drk_dc(num_rk+i,this%indices_rk%cols(num_rk+this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%col_1)=drk_dc_loc !> Commented: would store gradient in Jacobian (workaround: chapuza)
        !deallocate(drk_dc_loc) !> Commented: would free gradient array
    end do !> End mineral kinetic reactions loop
    num_rk=num_rk+this%solid_chemistry%mineral_zone%num_minerals_kin !> Update total reaction counter: add number of kinetic minerals to num_rk [-]
end subroutine !> End of compute_rk_new subroutine