!> @file aq_phase_m.f90
!> @brief Aqueous phase module for reactive transport modeling
!> @details This module defines the aqueous phase class (subclass of phase_c) that:
!> - Contains all aqueous species in a solution
!> - Computes log₁₀ activity coefficients using various models (Debye-Hückel, Davies, etc.)
!> - Computes Jacobian of log₁₀ activity coefficients with respect to log₁₀ concentrations
!> - Rearranges aqueous species into variable and constant activity groups
!> - Tracks important species indices (water, H⁺, OH⁻, HCO₃⁻, CO₃²⁻)
!>
!> @par Activity Coefficients:
!> Activity coefficients γᵢ relate activity aᵢ to concentration cᵢ: aᵢ = γᵢ·cᵢ
!> For aqueous ions, log₁₀(γᵢ) depends on ionic strength I according to models like:
!> - Debye-Hückel: log₁₀(γᵢ) = -A·zᵢ²·√I / (1 + a₀·B·√I)
!> - Davies: log₁₀(γᵢ) = -A·zᵢ²·[√I/(1+√I) - 0.3·I]
!>
!> @par Jacobian Computation:
!> For speciation calculations, we need ∂log₁₀(γᵢ)/∂log₁₀(cⱼ) which involves:
!> dlog₁₀(γᵢ)/dlog₁₀(cⱼ) = (dlog₁₀(γᵢ)/dI) · (dI/dcⱼ) · (dcⱼ/dlog₁₀(cⱼ))
!>
!> @see phase_m For base phase class
!> @see aq_species_m For aqueous species class
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module aq_phase_m
    use phase_m, only : phase_c !< Import base phase class for inheritance
    use reaction_m, only : reaction_c !< Import reaction class for reaction indexing
    use aq_species_m, only : aq_species_c !< Import aqueous species classes
    use species_m, only : species_c !< Import base species class
    use arrays_m, only : diag_matrix_c !< Import diagonal matrix class for Jacobian calculations
    use params_aq_sol_m, only : params_aq_sol_s !< Import aqueous solution parameters (Debye-Hückel constants)
    !use vectors_m !< Reserved for future vector operations
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Default visibility is private; expose only intended interfaces
    
    !> @brief Aqueous phase class extending base phase
    !> @details This type represents an aqueous solution containing dissolved species.
    !> It manages species organization, activity coefficient calculations, and important species tracking.
    !>
    !> @par Species Organization:
    !> Species are arranged by activity behavior:
    !> 1. Variable activity species (most ions and complexes)
    !> 2. Constant activity species (e.g., H₂O in dilute solutions)
    !>
    !> @par Important Species Tracking:
    !> The class tracks indices of key species for pH and carbonate chemistry:
    !> - Water (H₂O): Solvent, usually constant activity
    !> - H⁺: Proton, controls pH
    !> - OH⁻: Hydroxide, water dissociation product
    !> - HCO₃⁻: Bicarbonate, carbonate system
    !> - CO₃²⁻: Carbonate, carbonate system
    !>
    !> @extends phase_c
    type, public, extends(phase_c) :: aq_phase_c !< Aqueous phase class
        integer(kind=4) :: num_aq_complexes=0 !< Number of aqueous complexes (secondary species) [-] (default: 0)
        type(aq_species_c), allocatable :: aq_species(:) !< Array of aqueous species in solution (size: num_species) [-]
        integer(kind=4) :: wat_flag=0 !< ideal water presence flag: 1 if ideal water is present as a species, 0 otherwise [-] (default: 0)
        integer(kind=4) :: ind_wat=0 !< Index of water (H₂O) in aq_species array [-] (0 if absent)
        integer(kind=4), allocatable :: ind_diss_solids(:) !< Indices of dissolved solids (non-water species) in aq_species [-]
        integer(kind=4) :: ind_prot=0 !< Index of proton (H⁺) in aq_species array [-] (0 if absent)
        integer(kind=4) :: ind_oh=0 !< Index of hydroxide (OH⁻) in aq_species array [-] (0 if absent)
        integer(kind=4) :: ind_bicarb=0 !< Index of bicarbonate (HCO₃⁻) in aq_species array [-] (0 if absent)
        integer(kind=4) :: ind_carb=0 !< Index of carbonate (CO₃²⁻) in aq_species array [-] (0 if absent)
    contains
    !> @name Setter Methods
    !> @{
        procedure :: set_ind_prot !< Set index of proton (H⁺) in species array
        procedure :: set_ind_wat !< Set index of water (H₂O) in species array
        procedure :: set_ind_diss_solids !< Set indices of dissolved solids (non-water species)
        procedure :: set_indices_aq_phase !< Automatically detect and set indices of important species (H₂O, H⁺, OH⁻, HCO₃⁻, CO₃²⁻)
        procedure :: set_num_aq_complexes !< Set number of aqueous complexes (secondary species)
        procedure :: set_aq_species !< Set entire array of aqueous species
        procedure :: set_single_aq_species !< Set a single aqueous species at specified index
    !> @}
    !> @name Allocation Methods
    !> @{
        procedure :: allocate_aq_species !< Allocate aqueous species array
        procedure :: allocate_ind_diss_solids !< Allocate dissolved solids indices array
    !> @}
    !> @name Computation Methods
    !> @{
        procedure :: compute_log_act_coeffs_aq_phase !< Compute log₁₀ activity coefficients for all variable activity species
        procedure :: compute_log_Jacobian_act_coeffs_aq_phase !< Compute Jacobian ∂log₁₀(γ)/∂log₁₀(c) for speciation calculations
    !> @}
    !> @name Reorganization Methods
    !> @{
        procedure :: rearrange_aq_species !< Rearrange species: variable activity first, then constant activity
    !> @}
    !> @name Query Methods
    !> @{
        procedure :: is_species_in_aq_phase !< Check if a species is present in aqueous phase and return its index
        procedure :: is_water_in_aq_phase !< Check if water (H₂O) is present in aqueous phase
    !> @}
    !> @name Copy Methods
    !> @{
        procedure :: copy_attributes !< Copy all attributes from another aqueous phase object
    !> @}
    !> @name Getter Methods
    !> @{
        procedure :: get_valences !< Get array of valences (charges) for all aqueous species
        procedure :: get_indices_reaction   !< Get indices of species participating in specific reactions
        procedure :: get_aq_species_index_by_name !< Get index of aqueous species by name
    !> @}
    end type
    
    !> @brief Reserved interface block for future extensions
    !> @details Currently contains commented-out interface declarations.
    !> May be used for operator overloading or generic interfaces in future.
    interface
    
        function get_indices_reaction(this,reaction) result(indices)
            import aq_phase_c
            import reaction_c
            implicit none
            class(aq_phase_c), intent(in) :: this
            class(reaction_c), intent(in) :: reaction
            integer(kind=4), allocatable :: indices(:)
        end function
        ! subroutine compute_log_act_coeffs_aq_phase(this,ionic_strength,params_aq_sol,log_act_coeffs)
        !     import aq_phase_c
        !     import params_aq_sol_s
        !     implicit none
        !     class(aq_phase_c) :: this
        !     real(kind=8), intent(in) :: ionic_strength
        !     class(params_aq_sol_s), intent(in) :: params_aq_sol
        !     real(kind=8), intent(out) :: log_act_coeffs(:) !> must be allocated
        ! end subroutine
        
        ! subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,out_prod,conc,log_Jacobian_act_coeffs)
        !     import aq_phase_c
        !     implicit none
        !     class(aq_phase_c) :: this
        !     real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_d_I and z^2
        !     real(kind=8), intent(in) :: conc(:) !> concentration of species in a given target
        !     real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> must be allocated
        ! end subroutine        
        
    end interface
    
    contains
        
        !> @brief Set number of aqueous complexes (secondary species)
        !> @details Assigns the count of aqueous complexes formed by association reactions.
        !> Complexes are secondary species derived from primary species via equilibrium reactions.
        !>
        !> @par Example:
        !> If the system has primary species [Ca²⁺, CO₃²⁻, H⁺] and forms complex [CaCO₃⁰],
        !> then num_aq_complexes = 1.
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] num_aq_complexes Number of aqueous complexes [-]
        subroutine set_num_aq_complexes(this,num_aq_complexes)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
            integer(kind=4), intent(in) :: num_aq_complexes !< Number of aqueous complexes (secondary species) [-]
            this%num_aq_complexes=num_aq_complexes !< Assign number of complexes [-]
        end subroutine
        
        !> @brief Set index of proton (H⁺) in aqueous species array
        !> @details Assigns the array position of H⁺ for pH calculations and acid-base chemistry.
        !> The proton concentration controls pH via: pH = -log₁₀([H⁺])
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] ind_proton Index of H⁺ in aq_species array [-]
        subroutine set_ind_prot(this,ind_proton)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           integer(kind=4), intent(in) :: ind_proton !< Index of proton in aq_species array [-]
           this%ind_prot=ind_proton !< Assign proton index [-]
        end subroutine
        
        !> @brief Set index of water (H₂O) in aqueous species array
        !> @details Assigns the array position of water, which is typically treated as constant
        !> activity (aH₂O = 1) in dilute aqueous solutions.
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] ind_wat Index of H₂O in aq_species array [-]
        subroutine set_ind_wat(this,ind_wat)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           integer(kind=4), intent(in) :: ind_wat !< Index of water in aq_species array [-]
           this%ind_wat=ind_wat !< Assign water index [-]
        end subroutine
        
        !> @brief Allocate aqueous species array
        !> @details Allocates memory for the aq_species array. If already allocated, deallocates first.
        !> The array size is set from the optional argument or uses the existing num_species value.
        !>
        !> @par Usage:
        !> - With argument: call aq_phase%allocate_aq_species(10)  ! Allocate for 10 species
        !> - Without: call aq_phase%allocate_aq_species()  ! Use existing num_species
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] num_species (Optional) Number of aqueous species to allocate [-]
        subroutine allocate_aq_species(this,num_species)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
            integer(kind=4), intent(in), optional :: num_species !< Number of species (optional) [-]
            
            !> If num_species provided, update the count
            if (present(num_species)) then
                this%num_species=num_species !< Set number of species [-]
            end if
            !> Deallocate if already allocated (reallocate with new size)
            if (allocated(this%aq_species)) then
                deallocate(this%aq_species) !< Free existing memory
            end if
            !> Allocate array for aqueous species
            allocate(this%aq_species(this%num_species)) !< Allocate array of size num_species
        end subroutine
        
        !> @brief Allocate dissolved solids indices array
        !> @details Allocates memory for indices of dissolved species (excluding water).
        !> Array size = num_species - wat_flag (where wat_flag = 1 if water present, 0 otherwise).
        !>
        !> @par Example:
        !> If system has 10 species including water: size = 10 - 1 = 9 dissolved solids
        !>
        !> @param[inout] this Aqueous phase object
        subroutine allocate_ind_diss_solids(this)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           !> Allocate indices array for non-water species only
           allocate(this%ind_diss_solids(this%num_species-this%wat_flag)) !< Size excludes water if present
        end subroutine

        
        !> @brief Set entire array of aqueous species
        !> @details Assigns an array of aqueous species to the phase and automatically
        !> updates num_species based on the array size.
        !>
        !> @par Note:
        !> This will overwrite any existing species in the phase.
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] aq_species Array of aqueous species to assign [-]
        subroutine set_aq_species(this,aq_species)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
            class(aq_species_c), intent(in) :: aq_species(:) !< Array of aqueous species [-]
            this%aq_species=aq_species !< Assign species array [-]
            this%num_species=size(aq_species) !< Update species count from array size [-]
        end subroutine
        
        !> @brief Set a single aqueous species at specified index
        !> @details Assigns one species to a specific position in the aq_species array.
        !> Performs validation checks to ensure array is allocated and index is valid.
        !>
        !> @par Validation:
        !> - Array must be allocated before calling
        !> - Index must be positive (>= 1)
        !> - Index must not exceed num_species
        !>
        !> @param[inout] this Aqueous phase object
        !> @param[in] aq_species Single aqueous species to assign [-]
        !> @param[in] index Array position where species should be placed (1-indexed) [-]
        subroutine set_single_aq_species(this,aq_species,index)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
            class(aq_species_c), intent(in) :: aq_species !< Single aqueous species to assign [-]
            integer(kind=4), intent(in) :: index !< Array index (1-based) [-]
            
            !> Validate that species array is allocated
            if (.not. allocated(this%aq_species)) then
                error stop "Aqueous species must be allocated" !< Terminate if array not allocated
            !> Validate index is not too large
            else if (index>this%num_species) then
                error stop "Index is higher than number of aqueous species" !< Terminate if index out of bounds
            !> Validate index is positive
            else if (index<1) then
                error stop "Index must be positive" !< Terminate if index < 1
            !> All checks passed, assign species
            else
                this%aq_species(index)=aq_species !< Assign species at specified index [-]
            end if
        end subroutine
        
       !> @brief Rearrange aqueous species by activity behavior
       !> @details Reorganizes the aq_species array to group species by activity type:
       !> - **Variable activity species** (ions, complexes): Placed at beginning of array
       !> - **Constant activity species** (e.g., H₂O in dilute solutions): Placed at end of array
       !>
       !> @par Algorithm:
       !> 1. Variable activity non-water species → indices 1, 2, 3, ...
       !> 2. Variable activity water → next available index (if present)
       !> 3. Constant activity species → indices from end backwards (..., n-2, n-1, n)
       !>
       !> @par Purpose:
       !> This arrangement optimizes speciation calculations by separating species that
       !> require activity coefficient updates (variable) from those that don't (constant).
       !>
       !> @param[inout] this Aqueous phase object
       subroutine rearrange_aq_species(this)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           
           type(aq_species_c), allocatable :: aux(:) !< Temporary array to store original species order [-]
           integer(kind=4) :: i !< Loop counter over all species [-]
           integer(kind=4) :: j !< Counter for variable activity species positions (forward from 1) [-]
           integer(kind=4) :: k !< Counter for constant activity species (backward from end) [-]
           integer(kind=4) :: l !< Counter for dissolved solids with variable activity [-]
           integer(kind=4) :: m !< Counter for dissolved solids with constant activity [-]
           
           aux=this%aq_species !< Save original species array (temporary copy) [-]
           j=1 !< Initialize forward counter: variable activity species start at index 1 [-]
           l=1 !< Initialize dissolved solids (variable) counter [-]
           k=0 !< Initialize backward counter: constant activity species from end [-]
           m=0 !< Initialize dissolved solids (constant) counter [-]
           
           !> Loop through all species in original order
           do i=1,this%num_species
               !> Case 1: Variable activity non-water species (most ions and complexes)\n
               if (aux(i)%cst_act_flag .eqv. .false. .and. aux(i)%name/='h2o(p)' .and. aux(i)%name/='h2o') then
                    this%ind_diss_solids(l)=j !< Record index in dissolved solids array [-]
                    l=l+1 !< Increment dissolved solids counter [-]
                    call this%aq_species(j)%copy_species(aux(i)) !< Place at forward position j [-]
                    if (j<this%num_species) then
                        j=j+1 !< Move forward counter to next position [-]
                    end if
               !> Case 2: Constant activity non-water species\n
               else if (aux(i)%name/='h2o(p)' .and. aux(i)%name/='h2o') then
                    this%ind_diss_solids(this%num_species-m)=this%num_species-k !< Record backward index [-]
                    m=m+1 !< Increment constant dissolved solids counter [-]
                    call this%aq_species(this%num_species-k)%copy_species(aux(i)) !< Place at backward position [-]
                    k=k+1 !< Move backward counter [-]
               !> Case 3: Water with variable activity (rare, occurs at high concentrations)\n
               else if (aux(i)%cst_act_flag .eqv. .false.) then
                    call this%set_ind_wat(j) !< Set water index to current forward position [-]
                    call this%aq_species(j)%copy_species(aux(i)) !< Place water at forward position [-]
                    if (j<this%num_species) then
                        j=j+1 !< Move forward counter [-]
                    end if
               !> Case 4: Water with constant activity (typical: aH₂O = 1 in dilute solutions)\n
               else
                    call this%set_ind_wat(this%num_species-k) !< Set water index to backward position [-]
                    call this%aq_species(this%num_species-k)%copy_species(aux(i)) !< Place water at backward position [-]
                    k=k+1 !< Move backward counter [-]
               end if
           end do
       end subroutine
        
        ! subroutine append_aq_species(this,aq_species)
        !     implicit none
        !     class(aq_phase_c) :: this
        !     type(aq_species_c), intent(in) :: aq_species
            
        !     integer(kind=4) :: i
        !     type(aq_species_c), allocatable :: aux_array(:)
            
        !     aux_array=this%aq_species
        !     deallocate(this%aq_species)
        !     if (size(aux_array)>0) then
        !         call this%allocate_aq_species(size(aux_array)+1)
        !         do i=1,this%num_species-1
        !             this%aq_species(i)=aux_array(i)
        !         end do
        !         this%aq_species(this%num_species)=aq_species
        !     else
        !         call this%allocate_aq_species(1)
        !         call this%set_aq_species([aq_species]) !> chapuza
        !     end if
        ! end subroutine
        
        !> @brief Check if a species is present in aqueous phase
        !> @details Searches for a species by name in the aq_species array.
        !> Returns a logical flag indicating presence and optionally returns the array index.
        !>
        !> @par Usage:
        !> ```fortran
        !> logical :: present
        !> integer :: idx
        !> call aq_phase%is_species_in_aq_phase(my_species, present, idx)
        !> if (present) then
        !>     ! Species found at position idx
        !> end if
        !> ```
        !>
        !> @param[in] this Aqueous phase object
        !> @param[in] species Species to search for (searches by name)
        !> @param[out] flag True if species found, False otherwise [-]
        !> @param[out] species_ind (Optional) Index of species in aq_species array (0 if not found) [-]
        subroutine is_species_in_aq_phase(this,species,flag,species_ind)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
            class(species_c), intent(in) :: species !< Species to search for [-]
            logical, intent(out) :: flag !< True if found, False otherwise [-]
            integer(kind=4), intent(out), optional :: species_ind !< Index if found (optional) [-]
            
            integer(kind=4) :: i !< Loop counter over species [-]
            
            flag=.false. !< Initialize: assume species not present [-]
            !> Initialize optional index to 0 (not found)
            if (present(species_ind)) then
                species_ind=0 !< Set to 0 if species not found [-]
            end if
            !> Search through all species in aqueous phase
            do i=1,this%num_species
                !> Compare species names (case-sensitive)
                if (species%name==this%aq_species(i)%name) then
                    flag=.true. !< Species found [-]
                    !if (species%name=='h2o') then
                    !    this%ind_wat=i  !< Commented: Could auto-set water index
                    !end if
                    !> Return index if requested
                    if (present(species_ind)) then
                        species_ind=i !< Return array position [-]
                    end if
                    exit !< Exit loop once found (names should be unique)
                end if
            end do
        end subroutine
        
        !> @brief Check if water (H₂O) is present in aqueous phase
        !> @details Convenience wrapper that searches specifically for water by name 'h2o'.
        !> Uses is_species_in_aq_phase internally.
        !>
        !> @par Note:
        !> Water is often present in aqueous systems as the solvent, typically with
        !> constant activity (aH₂O = 1) in dilute solutions.
        !>
        !> @param[in] this Aqueous phase object
        !> @param[out] flag True if water found, False otherwise [-]
        subroutine is_water_in_aq_phase(this,flag)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
            logical, intent(out) :: flag !< True if water present, False otherwise [-]
            
            type(aq_species_c) :: water !< Temporary water species for search [-]

            call water%set_name('h2o') !< Set name to search for (lowercase convention) [-]
            call this%is_species_in_aq_phase(water,flag) !< Search for water species [-]
            if (flag .eqv. .false.) then
               call water%set_name('h2o(p)') !< Set name to search for (lowercase convention) [-]
               call this%is_species_in_aq_phase(water,flag) !< Search for water species [-]
            end if
        end subroutine
        
        !> @brief Automatically detect and set indices of important species
        !> @details Scans the aq_species array and identifies key species by name,
        !> storing their indices for quick access in chemistry calculations:
        !> - H₂O (water): Solvent
        !> - H⁺ (proton): pH control
        !> - OH⁻ (hydroxide): Water dissociation product
        !> - HCO₃⁻ (bicarbonate): Carbonate system
        !> - CO₃²⁻ (carbonate): Carbonate system
        !>
        !> @par Name Matching:
        !> Uses lowercase names with charge notation:
        !> - 'h2o(p)' for water
        !> - 'h+' for proton
        !> - 'oh-' for hydroxide
        !> - 'hco3-' for bicarbonate
        !> - 'co3-2' for carbonate
        !>
        !> @param[inout] this Aqueous phase object
        subroutine set_indices_aq_phase(this)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           
           integer(kind=4) :: i !< Loop counter over species [-]
           
           !> Loop through all species to identify important ones by name
           do i=1,this%num_species
               !> Check for water (phase notation)
               if (this%aq_species(i)%name=='h2o(p)') then
                  this%ind_wat=i !< Store water index [-]
               !> Check for proton
               else if (this%aq_species(i)%name=='h+') then
                   this%ind_prot=i !< Store proton index [-]
               !> Check for hydroxide
               else if (this%aq_species(i)%name=='oh-') then
                   this%ind_oh=i !< Store hydroxide index [-]
               !> Check for bicarbonate
               else if (this%aq_species(i)%name=='hco3-') then
                   this%ind_bicarb=i !< Store bicarbonate index [-]
               !> Check for carbonate
               else if (this%aq_species(i)%name=='co3-2') then
                   this%ind_carb=i !< Store carbonate index [-]
               !> Species not in important list
               else
                   continue !< Skip to next species
               end if
           end do
        end subroutine
        
        !> @brief Set indices of dissolved solids (non-water species)
        !> @details Populates the ind_diss_solids array with indices of all species
        !> except water. This array is used to iterate over dissolved species only.
        !>
        !> @par Array Size:
        !> ind_diss_solids has size = num_species - wat_flag
        !> - If water present (wat_flag=1): size = num_species - 1
        !> - If no water (wat_flag=0): size = num_species
        !>
        !> @par Example:
        !> If aq_species = [Na⁺, Cl⁻, H₂O, Ca²⁺] then ind_diss_solids = [1, 2, 4]
        !>
        !> @param[inout] this Aqueous phase object
        subroutine set_ind_diss_solids(this)
           implicit none !< Enforce explicit variable declarations
           class(aq_phase_c), intent(inout) :: this !< Aqueous phase object [-]
           
           integer(kind=4) :: i !< Loop counter over all species [-]
           integer(kind=4) :: j !< Counter for dissolved solids array [-]
           
           j=1 !< Initialize dissolved solids counter [-]
           !> Loop through all species
           do i=1,this%num_species
               !> Check if species is not water
               if (this%aq_species(i)%name/='h2o') then
                   this%ind_diss_solids(j)=i !< Store index of non-water species [-]
                   !> Check if more dissolved solids remain
                   if (j<this%num_species-this%wat_flag) then
                       j=j+1 !< Increment dissolved solids counter [-]
                   else
                       exit !< All dissolved solids found, exit loop
                   end if
               end if
           end do
        end subroutine
                
        !> @brief Copy all attributes from another aqueous phase
        !> @details Deep copy of all aqueous phase attributes including:
        !> - All species (via copy_species)
        !> - Water flag and index
        !> - Number of complexes
        !> - Dissolved solids indices
        !>
        !> @par Allocation:
        !> If this%aq_species not yet allocated, allocates it to match source size.
        !> If already allocated, uses existing allocation.
        !>
        !> @param[inout] this Destination aqueous phase object
        !> @param[in] aq_phase Source aqueous phase to copy from
        subroutine copy_attributes(this,aq_phase)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(inout) :: this !< Destination aqueous phase [-]
            class(aq_phase_c), intent(in) :: aq_phase !< Source aqueous phase to copy [-]
            
            integer(KIND=4) :: i !< Loop counter over species [-]
            
            !> Check if species array needs allocation
            if (allocated(this%aq_species)) then
                continue !< Array already allocated, proceed with copy
            else
                call this%allocate_aq_species(aq_phase%num_species) !< Allocate to match source size
            end if
            !> Copy scalar attributes
            this%wat_flag=aq_phase%wat_flag !< Copy water presence flag [-]
            this%num_aq_complexes=aq_phase%num_aq_complexes !< Copy number of complexes [-]
            this%ind_wat=aq_phase%ind_wat !< Copy water index [-]
            !this%ind_proton=aq_phase%ind_proton !< Commented: proton index (not currently used)
            this%ind_diss_solids=aq_phase%ind_diss_solids !< Copy dissolved solids indices array [-]
            !> Deep copy each species individually
            do i=1,this%num_species
                call this%aq_species(i)%copy_species(aq_phase%aq_species(i)) !< Copy species i [-]
            end do
        end subroutine
        
        !> @brief Get array of valences (charges) for all aqueous species
        !> @details Extracts the valence (charge number) of each species and returns
        !> them as an integer array. Valences are used in activity coefficient calculations.
        !>
        !> @par Valence Examples:
        !> - Na⁺: valence = +1
        !> - Ca²⁺: valence = +2
        !> - Cl⁻: valence = -1
        !> - SO₄²⁻: valence = -2
        !> - H₂O: valence = 0 (neutral)
        !>
        !> @param[in] this Aqueous phase object
        !> @return valences Array of valences for all species (size: num_species) [-]
        function get_valences(this) result(valences)
            implicit none !< Enforce explicit variable declarations
            class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
            integer(kind=4), allocatable :: valences(:) !< Output array of valences [-]
            
            integer(kind=4) :: i !< Loop counter over species [-]
            
            !> Allocate result array for all species
            allocate(valences(this%num_species)) !< Allocate array of size num_species
            !> Extract valence from each species
            do i=1,this%num_species
                valences(i)=this%aq_species(i)%valence !< Copy valence of species i [-]
            end do
        end function

         !> @brief Compute log₁₀ activity coefficients for dissolved solids
         !> @details Calculates log₁₀(γᵢ) for all dissolved species (excluding water) using
         !> the generalized activity coefficient model:
         !> log₁₀(γᵢ) = -α·A·zᵢ²·√I / (1 + β·√I) + γ·I
         !> where α, β, γ are model-specific parameters stored in params_act_coeff.
         !>
         !> @par Supported Models (via α, β, γ):
         !> - Debye-Hückel restricted: α=1, β=0, γ=0
         !> - Debye-Hückel extended: α=1, β=a₀·B, γ=0
         !> - Davies: α=1, β=1, γ=0.3·A·z²
         !> - Truesdell-Jones: α=1, β=B·aTJ, γ=bTJ
         !> - b-dot (Helgeson): α=1, β=B·aTJ, γ=ḃ(T)
         !>
         !> @par Ionic Strength:
         !> I = 0.5·Σ(cᵢ·zᵢ²) where cᵢ is concentration [mol/L] and zᵢ is valence
         !>
         !> @par Water Treatment:
         !> Water is typically excluded (constant activity aH₂O = 1 in dilute solutions).
         !> Only dissolved solids (ions and complexes) require activity corrections.
         !>
         !> @param[in] this Aqueous phase object
         !> @param[in] ionic_strength Ionic strength I [mol/L]
         !> @param[in] params_aq_sol Aqueous solution parameters (Debye-Hückel constants A, B)
         !> @param[out] log_act_coeffs Array of log₁₀(γᵢ) values (must be pre-allocated, size: num_species) [-]
subroutine compute_log_act_coeffs_aq_phase(this,ionic_strength,params_aq_sol,log_act_coeffs)
    implicit none !< Enforce explicit variable declarations
    class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
    real(kind=8), intent(in) :: ionic_strength !< Ionic strength I [mol/L]
    class(params_aq_sol_s), intent(in) :: params_aq_sol !< Debye-Hückel parameters (A, B) [-]
    real(kind=8), intent(out) :: log_act_coeffs(:) !< Output: log₁₀ activity coefficients (pre-allocated) [-]
    
    integer(kind=4) :: i !< Loop counter over dissolved solids [-]
    
    !> Commented: Model selection based on ionic strength (adaptive model)
    !if (ionic_strength<1d-2) then
    !    p_compute_log_act_coeff=>Debye_Huckel_restr  !< Low I: use restricted DH
    !else if (ionic_strength<=1d-1) then
    !    p_compute_log_act_coeff=>Debye_Huckel_ampl   !< Medium I: use extended DH
    !else if (ionic_strength<=7d-1) then
    !    p_compute_log_act_coeff=>Davies              !< High I: use Davies
    !else
    !    error stop                                    !< Very high I: not supported
    !end if
    !p_compute_log_act_coeff=>Davies !< Commented: Force Davies model (temporary)
    
    !> Loop over all dissolved solids (non-water species)
    do i=1,this%num_species-this%wat_flag
        !call p_compute_log_act_coeff(this%aq_species(i),ionic_strength,log_act_coeffs(i))  !< Commented: function pointer approach
        
        !> Compute log₁₀(γᵢ) using generalized formula:
        !> log₁₀(γᵢ) = -α·A·zᵢ²·√I / (1 + β·√I) + γ·I
        log_act_coeffs(this%ind_diss_solids(i))=-this%aq_species(this%ind_diss_solids(i))%params_act_coeff%alpha*(params_aq_sol%A*&
        this%aq_species(this%ind_diss_solids(i))%valence**2)*sqrt(ionic_strength)/(1d0+this%aq_species(&
        this%ind_diss_solids(i))%params_act_coeff%beta*sqrt(ionic_strength)) +&
         this%aq_species(this%ind_diss_solids(i))%params_act_coeff%gamma*ionic_strength
         !< First term: -α·A·z²·√I/(1+β·√I)  [Debye-Hückel with extensions]
         !< Second term: +γ·I  [Linear ionic strength correction]
    end do
    !log_act_coeffs(this%ind_wat)  !< Commented: Water activity coefficient (typically log₁₀(γH₂O) = 0)
end subroutine

!> @brief Compute Jacobian of log₁₀ activity coefficients with respect to log₁₀ concentrations
!> @details Calculates the matrix ∂log₁₀(γᵢ)/∂log₁₀(cⱼ) required for Newton-Raphson
!> speciation calculations. This Jacobian accounts for how changing concentration of
!> species j affects the activity coefficient of species i through ionic strength.
!>
!> @par Mathematical Derivation:
!> The Jacobian is computed via chain rule:
!> ∂log₁₀(γᵢ)/∂log₁₀(cⱼ) = (∂log₁₀(γᵢ)/∂I) · (∂I/∂cⱼ) · (∂cⱼ/∂log₁₀(cⱼ))
!>
!> Since I = 0.5·Σ(cₖ·zₖ²):
!> ∂I/∂cⱼ = 0.5·zⱼ²
!>
!> And ∂cⱼ/∂log₁₀(cⱼ) = ln(10)·cⱼ:
!> ∂log₁₀(γᵢ)/∂log₁₀(cⱼ) = (∂log₁₀(γᵢ)/∂I) · (0.5·zⱼ²) · (ln(10)·cⱼ)
!>
!> @par Matrix Form:
!> log_Jacobian_act_coeffs = 0.5·ln(10)·diag(c)·out_prod
!> where out_prod = [∂log₁₀(γᵢ)/∂I]·[zⱼ²]ᵀ (outer product)
!>
!> @par Assumption:
!> All primary species are aqueous (valid for typical speciation problems)
!>
!> @param[in] this Aqueous phase object
!> @param[in] out_prod Outer product matrix: [∂log₁₀(γ)/∂I] ⊗ [z²] (size: n×n) [-]
!> @param[in] conc Subset of aqueous species concentrations [mol/L]
!> @param[out] log_Jacobian_act_coeffs Jacobian matrix ∂log₁₀(γ)/∂log₁₀(c) (must be pre-allocated) [-]
subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,out_prod,conc,log_Jacobian_act_coeffs)
    implicit none !< Enforce explicit variable declarations
    class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
    real(kind=8), intent(in) :: out_prod(:,:) !< Outer product: [∂log₁₀(γ)/∂I] ⊗ [z²] [-]
    real(kind=8), intent(in) :: conc(:) !< Subset of aqueous species concentrations [mol/L]
    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !< Jacobian: ∂log₁₀(γ)/∂log₁₀(c) (pre-allocated) [-]

    type(diag_matrix_c) :: conc_diag !< Diagonal matrix with concentrations on diagonal [mol/L]
    
    !> Create diagonal matrix from concentration vector
    call conc_diag%set_diag_matrix(conc) !< diag(c) where c = [c₁, c₂, ..., cₙ] [mol/L]
    
    !> Compute Jacobian: J = 0.5·ln(10)·diag(c)·outer_product
    !> Factor breakdown:
    !> - 5d-1 = 0.5 (from ∂I/∂c = 0.5·z²)
    !> - log(1d1) = ln(10) (from ∂c/∂log₁₀(c) = ln(10)·c)
    !> - conc_diag%prod_mat_diag_mat(out_prod) = diag(c)·outer_product
    log_Jacobian_act_coeffs=5d-1*log(1d1)*conc_diag%prod_mat_diag_mat(out_prod)
     
end subroutine

subroutine get_aq_species_index_by_name(this,name,index)
    implicit none !< Enforce explicit variable declarations
    class(aq_phase_c), intent(in) :: this !< Aqueous phase object [-]
    character(len=*), intent(in) :: name !< Name of the aqueous species [-]
    integer(kind=4), intent(out) :: index !< Index of the aqueous species [-]
    
    integer(kind=4) :: i  !< Loop counter over species [-]
    !> Search for the species by name
    index = 0 !< Initialize index to 0 (not found)
    do i=1,this%num_species
        if (trim(this%aq_species(i)%name) == trim(name)) then
            index = i !< Species found, set index
            exit
        end if
    end do
    if (index == 0) then
        error stop "Aqueous species not found"
    end if
end subroutine

        
end module