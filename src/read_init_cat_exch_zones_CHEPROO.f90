!> @file read_init_cat_exch_zones_CHEPROO.f90
!> @brief Reads initial cation exchange zones from CHEPROO format input file
!> 
!> @details
!> This subroutine parses initial cation exchange zone definitions from a CHEPROO-formatted
!> geochemical input file. Cation exchange zones represent surface adsorption sites where
!> cations can exchange according to equilibrium relationships (e.g., Gaines-Thomas convention).
!>
!> **Algorithm Overview:**
!>
!> The subroutine operates in several phases:
!>
!> **Phase 1: Initialization**
!> - Read number of surface adsorption zones
!> - Allocate initial cation exchange zones array
!> - Set default Gaines-Thomas convention for exchange reactions
!>
!> **Phase 2: Read Each Cation Exchange Zone**
!> - Read zone index
!> - Read surface properties:
!>   * c_int: Internal concentration
!>   * c_ext: External concentration
!>   * spec_sorb_surf: Specific sorption surface area
!>   * CEC: Cation Exchange Capacity (equivalents/mass)
!> - Set up reactive zone with chemical system
!> - Configure cation exchange zone properties
!> - Allocate concentration, activity coefficient, and activity arrays
!> - Compute log activity coefficients for adsorbed cations using CEC
!>
!> **Phase 3: Integrate with Existing Reactive Zones**
!> - If gas reactive zones exist:
!>   * Create combined reactive zones (gas + cation exchange + gas+cat_exch)
!>   * Number of new zones = num_gas + num_cat_exch + (num_gas * num_cat_exch)
!> - If no gas zones:
!>   * Create simple cation exchange reactive zones
!>
!> **Cation Exchange Capacity (CEC):**
!> CEC represents the total number of exchangeable cation sites per unit mass of solid.
!> It is typically expressed in equivalents/kg (eq/kg) or meq/100g.
!>
!> **Gaines-Thomas Convention:**
!> Default convention for calculating selectivity coefficients in cation exchange reactions.
!> Activity coefficients depend on valences and equivalent fractions of adsorbed cations.
!>
!> @param[in,out] this Chemistry object containing all chemical system data
!> @param[in] unit File unit number for reading CHEPROO input file
!> @param[out] ndrz Number of adsorption reactive zones (output)
!>
!> @note The term "chapuza" (Spanish for "workaround/hack") appears in comments,
!>       indicating areas where the implementation uses temporary or non-ideal solutions
!>
!> @note Surface complexes are assumed to be in the same order as in the chemical system
!>
!> @warning Activity coefficients for adsorbed cations are computed using valences from
!>          the aqueous phase, which must be properly initialized before calling this routine
!>
!> @see chemistry_c For the main chemistry object structure
!> @see Gaines_Thomas_c For the Gaines-Thomas exchange convention
!> @see reactive_zone_c For reactive zone structure
!>
!> @author jordi
!> @date November 2025
!>
subroutine read_init_cat_exch_zones_CHEPROO(this,unit,ndrz)
!> ================================================================
!> MODULE IMPORTS
!> ================================================================
    use chemistry_m, only: chemistry_c
    !> chemistry classes: main chemistry, solid chemistry, reactive and mineral zones
    use Gaines_Thomas_m, only: Gaines_Thomas_c
    !> Gaines-Thomas convention for cation exchange
    use reactive_zone_m, only: reactive_zone_c !> reactive zone class
    use mineral_zone_m, only: mineral_zone_c !> mineral zone class
    implicit none !> no implicit variable typing
    
!> ================================================================
!> SUBROUTINE PARAMETERS
!> ================================================================
    class(chemistry_c) :: this !> chemistry object (contains all chemical system data)
    integer(kind=4), intent(in) :: unit !> file unit number for reading input file
    integer(kind=4), intent(out):: ndrz !> number of adsorption reactive zones (output)
    !type(solid_chemistry_c), intent(out), allocatable :: this%init_cat_exch_zones(:) 
        !> (COMMENTED) initial cation exchange zones (now part of chemistry object)
    !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
        !> (COMMENTED) optional reactive zones array
    
!> ================================================================
!> LOCAL VARIABLES - Counters and Indices
!> ================================================================
    integer(kind=4) :: i !> loop counter for reactive zones
    integer(kind=4) :: j !> loop counter for cation exchange zones
    integer(kind=4) :: k !> general loop counter (unused)
    integer(kind=4) :: ndtype !> total number of surface adsorption zone types
    integer(kind=4) :: idtype !> current adsorption zone type index
    integer(kind=4) :: num_min_rz !> number of existing (mineral) reactive zones read before the surface zones
    integer(kind=4) :: icon !> constraint icon (unused)
    integer(kind=4) :: num_ads_zones !> counter for adsorption zones processed
    integer(kind=4) :: num_mins_glob !> number of global minerals (unused)
    integer(kind=4) :: num_rz !> total number of reactive zones after combining
    integer(kind=4) :: min_idx !> index of the mineral zone merged into a cat-exch zone
    integer(kind=4) :: ns_merge !> number of surface complexes (merge helper)
    integer(kind=4) :: nmk_merge !> number of kinetic minerals (merge helper)
    integer(kind=4) :: off_merge !> offset of surface complexes in the combined solid arrays
    real(kind=8), allocatable :: old_conc(:) !> saved surface-complex concentrations (merge helper)
    real(kind=8), allocatable :: old_act(:) !> saved surface-complex activities (merge helper)
    real(kind=8), allocatable :: old_lac(:) !> saved surface-complex log activity coefficients (merge helper)
    
!> ================================================================
!> LOCAL VARIABLES - Dynamic Arrays
!> ================================================================
    integer(kind=4), allocatable :: valences(:) !> valences of aqueous species
    
!> ================================================================
!> LOCAL VARIABLES - Strings
!> ================================================================
    character(len=256) :: str !> temporary string for reading property labels
    character(len=256) :: constrain !> constraint name (unused)
    character(len=256) :: label !> section label (unused)
    
!> ================================================================
!> LOCAL VARIABLES - Real Numbers (Surface Properties)
!> ================================================================
    real(kind=8) :: c_int !> internal concentration
    real(kind=8) :: c_ext !> external concentration
    real(kind=8) :: spec_sorb_surf !> specific sorption surface area
    real(kind=8) :: CEC !> Cation Exchange Capacity (equivalents/mass)
    
!> ================================================================
!> LOCAL VARIABLES - Chemistry Objects
!> ================================================================
    type(reactive_zone_c) :: react_zone !> temporary reactive zone object
    type(mineral_zone_c) :: min_zone !> temporary mineral zone object (unused)
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> auxiliary array for copying existing reactive zones
    type(Gaines_Thomas_c), target :: Gaines_Thomas !> Gaines-Thomas convention object (target for pointer)

!> ================================================================
!> SECTION 1: Read number of adsorption zones and allocate arrays
!> ================================================================
    !> Read total number of surface adsorption zone types
    read(unit,*) ndtype !> number of surface adsorption zones to be defined

    !> Allocate initial cation exchange zones array in chemistry object
    call this%allocate_init_cat_exch_zones(ndtype) !> allocate array (size = ndtype)

    !> Set output parameter (number of adsorption reactive zones)
    ndrz=ndtype !> chapuza (copy ndtype to output parameter)
    
    ! call this%set_num_materials(ndtype) !> (COMMENTED) set number of materials
    ! call this%allocate_materials() !> (COMMENTED) allocate materials array

!> ================================================================
!> SECTION 2: Read and configure each cation exchange zone
!> ================================================================
!> Assumption: Surface complexes are in the same order as in the chemical system
!> Vamos a asumir que los complejos de superficie estan en el mismo orden que en el sistema quimico
    
    !> Initialize adsorption zones counter
    num_ads_zones=0 !> counter for processed adsorption zones
    
    do !> infinite loop - exit when all zones processed
        !> ----------------------------------------------------------------
        !> Read zone index
        !> ----------------------------------------------------------------
        read(unit,*) idtype !> read adsorption zone type index
        
        !> Validate zone index
        if (idtype<0 .or. idtype>ndtype) error stop !> index must be in range [0, ndtype]
        
        !> ----------------------------------------------------------------
        !> Read surface properties
        !> ----------------------------------------------------------------
        !> Read internal concentration
        read(unit,*) str, c_int !> read label and internal concentration value
        
        !> Read external concentration
        read(unit,*) str, c_ext !> read label and external concentration value
        
        !> Read specific sorption surface area
        read(unit,*) str, spec_sorb_surf !> read label and specific surface area
        
        !> Read Cation Exchange Capacity
        read(unit,*) str, CEC !> read label and CEC value (equivalents/mass)
        
        !> ----------------------------------------------------------------
        !> Set up reactive zone for this cation exchange zone
        !> ----------------------------------------------------------------
        !> Set chemical system in reactive zone
        call react_zone%set_chem_syst_react_zone(this%chem_syst)
            !> copy chemical system definition to reactive zone
        
        !> Configure as cation exchange zone
        call react_zone%set_cat_exch_zone()
            !> allocate and initialize cation exchange zone structure
        
        !> Set exchange convention to Gaines-Thomas
        react_zone%cat_exch_zone%convention=>Gaines_Thomas !> pointer assignment (by default)
        
        !> Set number of solid species
        call react_zone%set_num_solids()
            !> count solid species in cation exchange zone
        
        !> Set number of non-flowing species (= num_minerals + num_surf_compl
        !> + num_gases_eq). Without this call num_non_flow_species stays 0,
        !> which later forces set_speciation_alg_dimensions into its "no
        !> non-flow species" branch and yields num_eq_reactions = 0 for the
        !> water-type reactive zone (causing compute_inv_Se_2 to receive an
        !> empty / non-square submatrix).
        call react_zone%set_num_non_flow_species()
        
        !> Set non-flowing species
        call react_zone%set_ind_non_flow_species()
            !> identify species that don't flow (adsorbed species)
        
        !> Set control volume parameters
        call react_zone%set_CV_params(this%CV_params)
            !> copy control volume parameters from chemistry object
        
        !> ----------------------------------------------------------------
        !> Initialize cation exchange zone in chemistry object
        !> ----------------------------------------------------------------
        !> Set number of solids from surface complexes count
        this%init_cat_exch_zones(idtype)%num_solids=react_zone%cat_exch_zone%num_surf_compl
            !> number of solids = number of surface complexes
        
        !> Link reactive zone to cation exchange zone
        call this%init_cat_exch_zones(idtype)%set_reactive_zone(react_zone)
            !> copy reactive zone configuration
        
        !call min_zone%set_chem_syst_min_zone(this%chem_syst) 
            !> (COMMENTED) set chemical system for mineral zone by default
        
        !> Set dummy mineral zone (no minerals in cation exchange zones)
        call this%init_cat_exch_zones(idtype)%set_mineral_zone(this%min_zone_dummy)
            !> use empty mineral zone (cat exch zones don't have minerals)
        
        !> ----------------------------------------------------------------
        !> Allocate concentration and activity arrays
        !> ----------------------------------------------------------------
        !> Allocate solid concentrations array
        call this%init_cat_exch_zones(idtype)%allocate_conc_solids()
            !> allocate conc array (size = num_solids = num_surf_compl)
        
        !> Set Cation Exchange Capacity
        call this%init_cat_exch_zones(idtype)%set_CEC(CEC)
            !> store CEC value (equivalents/mass)
        
        !> Allocate equivalents array
        call this%init_cat_exch_zones(idtype)%allocate_equivalents()
            !> allocate array for equivalent fractions
        
        !> Allocate log activity coefficients array
        call this%init_cat_exch_zones(idtype)%allocate_log_act_coeffs_solid_chem()
            !> allocate log(gamma) array for solid species
        
        !> ----------------------------------------------------------------
        !> Compute activity coefficients for adsorbed cations
        !> ----------------------------------------------------------------
        !> Get valences from aqueous phase
        valences=this%chem_syst%aq_phase%get_valences()
            !> retrieve valences of all aqueous species
        
        !> Compute log activity coefficients for adsorbed cations
        call this%init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%compute_log_act_coeffs_ads_cats(valences(&
            this%init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%exch_cat_indices),CEC,&
            this%init_cat_exch_zones(idtype)%log_act_coeffs(&
            2:this%init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%num_surf_compl))
            !> compute log(gamma) using Gaines-Thomas convention
            !> inputs: valences of exchange cations, CEC
            !> output: log_act_coeffs array (indices 2 to num_surf_compl, index 1 is free site)
        
        !> Allocate activities array
        call this%init_cat_exch_zones(idtype)%allocate_activities()
            !> allocate activities array (a = gamma * c)
        
        !call this%init_cat_exch_zones(idtype)%set_conc_free_site() 
            !> (COMMENTED) set concentration of free site 'x-'
        !call this%init_cat_exch_zones(idtype)%compute_num_solids_solid_chem() 
            !> (COMMENTED) compute number of solid species
        !this%init_cat_exch_zones(idtype)%activities(1)=1d0-SUM(this%init_cat_exch_zones(idtype)%activities(2:this%init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%num_surf_compl))
            !> (COMMENTED) activity of free site = 1 - sum of adsorbed cation activities
        
        !> ----------------------------------------------------------------
        !> Increment counter and check for completion
        !> ----------------------------------------------------------------
        num_ads_zones=num_ads_zones+1 !> increment processed zones counter
        
        if (num_ads_zones==ndtype) exit !> if all zones processed, exit loop
    end do !> end cation exchange zones reading loop
    
!> ================================================================
!> SECTION 3: Integrate cation exchange zones with reactive zones
!> ================================================================
    !if (present(reactive_zones)) then !> (COMMENTED) if optional reactive_zones parameter provided
        
        if (allocated(this%reactive_zones)) then !> mineral reactive zones already exist (read before the surface zones)
            !> ----------------------------------------------------------------
            !> Case 1: keep the mineral-only zones and append combined mineral+exchange zones
            !> ----------------------------------------------------------------
            !> Pre-existing reactive zones are the MINERAL zones (gas zones are read AFTER the
            !> surface adsorption zones, so none exist yet). Keep them unchanged as the
            !> mineral-only solid zones (indices 1..num_min_rz) and APPEND one combined
            !> mineral+exchange reactive zone per surface adsorption zone (indices
            !> num_min_rz+1..num_min_rz+ndtype), so <root>_tar_sol.dat can select either the
            !> mineral-only zone (index<=num_min_rz) or the combined calcite+exchange zone
            !> (index>num_min_rz) via its "solid zone index".
            num_min_rz=size(this%reactive_zones)
            allocate(aux_react_zones(num_min_rz))
            do i=1,num_min_rz
                call aux_react_zones(i)%copy_react_zone(this%reactive_zones(i))
            end do
            num_rz=num_min_rz+ndtype
            call this%allocate_reactive_zones(num_rz)
            !> Mineral-only zones (unchanged)
            do i=1,num_min_rz
                call this%reactive_zones(i)%copy_react_zone(aux_react_zones(i))
            end do
            !> Combined mineral+exchange zones: base each on its matching mineral reactive zone
            !> (so it CONTAINS the equilibrium minerals), then add the surface complexes and
            !> recompute the speciation algebra so num_eq_reactions / the stoichiometric matrix
            !> cover BOTH the mineral equilibria and the exchange reactions. flag_comp=.true.
            !> matches the chemical system (constant-activity water species ordered last).
            do i=1,ndtype
                min_idx=min(i,num_min_rz) !> matching mineral zone
                if (aux_react_zones(min_idx)%num_minerals>0) then
                    call this%reactive_zones(num_min_rz+i)%copy_react_zone(aux_react_zones(min_idx))
                    call this%reactive_zones(num_min_rz+i)%set_cat_exch_zone(&
                        this%init_cat_exch_zones(i)%reactive_zone%cat_exch_zone)
                    call this%reactive_zones(num_min_rz+i)%set_num_solids()
                    call this%reactive_zones(num_min_rz+i)%set_num_non_flow_species()
                    call this%reactive_zones(num_min_rz+i)%set_ind_non_flow_species()
                    call this%reactive_zones(num_min_rz+i)%set_speciation_alg_dimensions(.true.)
                    call this%reactive_zones(num_min_rz+i)%set_ind_eq_reacts()
                    call this%reactive_zones(num_min_rz+i)%allocate_ind_var_act_species()
                    call this%reactive_zones(num_min_rz+i)%set_stoich_mat_react_zone()
                    call this%reactive_zones(num_min_rz+i)%set_ind_mins_stoich_mat()
                    call this%reactive_zones(num_min_rz+i)%set_ind_gases_stoich_mat()
                else
                    call this%reactive_zones(num_min_rz+i)%copy_react_zone(this%init_cat_exch_zones(i)%reactive_zone)
                end if
                !> Back-link the exchange solid chemistry to its combined reactive zone.
                call this%init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(num_min_rz+i))
            end do

        else !> no existing reactive zones
            !> ----------------------------------------------------------------
            !> Case 2: Create simple cation exchange reactive zones
            !> ----------------------------------------------------------------
            
            !allocate(this%reactive_zones(ndtype)) !> (COMMENTED) allocate array
            
            !> Allocate reactive zones array (one per cation exchange zone)
            call this%allocate_reactive_zones(ndtype) !> allocate array (size = ndtype)
            
            do i=1,ndtype !> loop over all cation exchange zones
                !> Copy cation exchange reactive zone
                call this%reactive_zones(i)%copy_react_zone(this%init_cat_exch_zones(i)%reactive_zone)
                    !> copy reactive zone from initial cat exch zone
                
                !> Back-link reactive zone to initial cat exch zone
                call this%init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(i))
                    !> chapuza (create bidirectional link)
            end do
        end if
    !end if !> (COMMENTED) end if present(reactive_zones)
    
!> ================================================================
!> SECTION 4: Build the combined solid chemistry (minerals + surface complexes)
!> ================================================================
!> The surface-adsorption reactive zone now CONTAINS the equilibrium minerals (merged in
!> SECTION 3). Here the matching solid chemistry (init_cat_exch_zones) is rebuilt so its
!> solid-species arrays follow the combined [minerals, surface complexes] layout used
!> throughout the code: minerals occupy the leading block (equilibrium then kinetic) and
!> the surface complexes the trailing num_surf_compl slots. Without this a target solid
!> that references a cation-exchange solid zone loses its mineral volumetric fractions and
!> kinetic mineral reactions (the output then shows all-zero vol_fract for those waters,
!> inconsistent with <root>_quim_loc.dat). Mineral zones are read before this routine, so
!> the mineral solid chemistries already live in this%materials(1..num_mineral_zones).
    if (this%num_mineral_zones>0) then
        do i=1,ndtype
            min_idx=min(i,this%num_mineral_zones) !> map cat-exch zone i to a mineral zone
            associate(sc=>this%init_cat_exch_zones(i), mm=>this%materials(min_idx))
            ns_merge=sc%reactive_zone%cat_exch_zone%num_surf_compl !> surface complexes
            !> Total minerals in the combined solid chemistry: equilibrium minerals (now part
            !> of the merged reactive zone) plus kinetic minerals (carried by the mineral zone).
            nmk_merge=sc%reactive_zone%num_minerals+mm%mineral_zone%num_minerals_kin
            if (nmk_merge>0) then
                !> Preserve the surface-complex sub-arrays (currently occupying indices 1:ns)
                old_conc=sc%concentrations(1:ns_merge)
                old_act=sc%activities(1:ns_merge)
                old_lac=sc%log_act_coeffs(1:ns_merge)
                !> Associate the mineral zone (holds the kinetic minerals) and recompute the
                !> total solid count = reactive-zone solids (eq minerals + surf compl) + kinetic minerals
                sc%mineral_zone=>mm%mineral_zone
                sc%num_solids=sc%reactive_zone%num_solids+mm%mineral_zone%num_minerals_kin
                !> Reallocate the solid-species arrays for the combined layout. The allocate
                !> routines only free their primary array, so drop the auxiliary buffers first
                !> to avoid "already allocated" runtime errors.
                if (allocated(sc%conc_old)) deallocate(sc%conc_old)
                if (allocated(sc%conc_old_old)) deallocate(sc%conc_old_old)
                if (allocated(sc%log_Jacobian_act_coeffs)) deallocate(sc%log_Jacobian_act_coeffs)
                call sc%allocate_conc_solids()
                call sc%allocate_activities()
                call sc%allocate_log_act_coeffs_solid_chem()
                !> Minerals occupy the leading block (pure phases: c=a=1, log gamma=0)
                sc%concentrations(1:nmk_merge)=1d0
                sc%activities(1:nmk_merge)=1d0
                sc%log_act_coeffs(1:nmk_merge)=0d0
                !> Surface complexes are stored at the end of the solid arrays
                off_merge=sc%num_solids-ns_merge
                sc%concentrations(off_merge+1:off_merge+ns_merge)=old_conc
                sc%activities(off_merge+1:off_merge+ns_merge)=old_act
                sc%log_act_coeffs(off_merge+1:off_merge+ns_merge)=old_lac
                !> Mineral volumetric fractions and reactive surfaces from the mineral zone
                sc%min_pore_vol_ratios=mm%min_pore_vol_ratios
                sc%react_surfaces=mm%react_surfaces
                !> Mineral variable/constant activity classification indices
                if (allocated(mm%var_act_species_indices)) &
                    sc%var_act_species_indices=mm%var_act_species_indices
                if (allocated(mm%cst_act_species_indices)) &
                    sc%cst_act_species_indices=mm%cst_act_species_indices
            end if
            end associate
        end do
    end if

!> ================================================================
!> SECTION 5: Post-processing (COMMENTED OUT)
!> ================================================================
!> Chapuza: Set materials in chemistry object
    !this%materials=this%init_cat_exch_zones
        !> (COMMENTED) assign initial cat exch zones to materials array
        
end subroutine