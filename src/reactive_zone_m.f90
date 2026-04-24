!> \file reactive_zone_m.f90
!> \brief Reactive zone module for spatial heterogeneity in geochemistry
!> \details
!>   This module defines reactive zones, which are abstract zones with:
!>   - Specific set of heterogeneous equilibrium reactions
!>   - Non-flowing (immobile) species (minerals, surface complexes, gases)
!>   
!>   Reactive zones are fundamental for:
!>   - Spatial variability in geochemical systems
!>   - Reactive transport modeling
!>   - Mineralogical heterogeneity
!>   - Biogeochemical zonation
!>   
!>   Key features:
!>   - Polymorphic species arrays (minerals, gases, surface complexes)
!>   - Stoichiometric matrices for mass action laws
!>   - Speciation algebra for equilibrium calculations
!>   - Cation exchange zones
!>   - Gas phases in equilibrium
!>   
!>   Non-flowing species organization:
!>   1. Constant activity minerals
!>   2. Constant activity gases
!>   3. Variable activity minerals
!>   4. Surface complexes (cation exchange)
!>   5. Variable activity gases
!>   
!>   Equilibrium reactions ordering:
!>   - If flag_comp = TRUE:
!>     * Constant activity minerals
!>     * Constant activity gases
!>     * Redox reactions
!>     * Aqueous complexation
!>     * Variable activity minerals
!>     * Cation exchange
!>     * Variable activity gases
!>   - If flag_comp = FALSE:
!>     * Redox reactions
!>     * Aqueous complexation
!>     * All minerals
!>     * Cation exchange
!>     * All gases
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module reactive_zone_m
    use chem_system_m, only: chem_system_c
    use species_m, only: species_c
    use surf_compl_m, only: cat_exch_zone_c
    use gas_phase_m, only: gas_phase_c
    use CV_params_m, only: CV_params_s
    use speciation_algebra_m, only: speciation_algebra_s
    use arrays_m, only: id_matrix
    implicit none
    save
    private !< Make all entities private by default
    public :: compare_react_zones !< Make the comparison function public
    !> \brief Reactive zone type for spatial geochemical heterogeneity
    !> \details
    !>   Defines a spatial region with specific geochemical characteristics.
    !>   
    !>   A reactive zone encapsulates:
    !>   - Set of equilibrium reactions active in the zone
    !>   - Non-flowing species (minerals, surface complexes, gases)
    !>   - Stoichiometric relationships
    !>   - Speciation algebra for equilibrium calculations
    !>   
    !>   Organization of non-flowing species:
    !>   1. Minerals with constant activity (infinite reservoir)
    !>   2. Gases with constant activity (buffered)
    !>   3. Minerals with variable activity (finite amount)
    !>   4. Surface complexes (cation exchange sites)
    !>   5. Gases with variable activity
    !>   
    !>   Applications:
    !>   - Reactive transport with mineralogical heterogeneity
    !>   - Biogeochemical zonation (redox zones)
    !>   - Spatial variability in weathering profiles
    !>   - Multiphase flow with chemistry
    !>   
    !>   Example usage:
    !>   ```fortran
    !>   type(reactive_zone_c) :: oxic_zone, anoxic_zone
    !>   ! Oxic zone: has O2, no sulfides
    !>   ! Anoxic zone: no O2, has sulfides
    !>   ```
    type, public :: reactive_zone_c
        integer(kind=4) :: id=0                                      !< Zone ID for labeling (non-negative)
        integer(kind=4) :: num_non_flow_species=0                !< Total number of immobile species
        integer(kind=4), allocatable :: ind_non_flow_species(:)      !< array of indices of non-flowing species in chemical system (same order as in chemical system)
        integer(kind=4) :: num_minerals=0                            !< Total number of minerals in equilibrium
        integer(kind=4) :: num_minerals_cst_act=0                    !< Number of constant activity minerals
        integer(kind=4) :: num_minerals_var_act=0                    !< Number of variable activity minerals
        integer(kind=4), allocatable :: ind_mins_chem_syst(:)        !< Mineral indices in chemical system (same order as in chemical system)
        type(cat_exch_zone_c) :: cat_exch_zone                            !< Cation exchange zone
        integer(kind=4) :: num_solids=0                              !< Number of solid phases (minerals + surf compl)
        type(gas_phase_c) :: gas_phase                               !< Gas phase in equilibrium
        !integer(kind=4) :: ind_gas_phase                           !< Number of equilibrium reactions in reactive zone
        !type(biofilm_c) :: biofilm                                   !< Biofilm phase
        real(kind=8), allocatable :: stoich_mat(:,:)                 !< Stoichiometric matrix for all reactions
        integer(kind=4), allocatable :: ind_mins_stoich_mat(:)       !< Mineral indices in stoich matrix (same order as ind_mins_chem_syst)
        integer(kind=4), allocatable :: ind_gases_stoich_mat(:)      !< Gas indices in stoich matrix
        integer(kind=4), allocatable :: ind_eq_reacts(:)             !< Heterogeneous equilibrium reaction indices in chemical system (same order as in chemical system)
        class(chem_system_c), pointer :: chem_syst                   !< Pointer to chemical system (shared)
        type(speciation_algebra_s) :: speciation_alg                 !< Speciation algebra object
        class(CV_params_s), pointer :: CV_params                     !< Convergence parameters for speciation/mixing
        real(kind=8), allocatable :: U_SkT_prod(:,:)                 !< Product U*S_k,nc^T for computations
        integer(kind=4), allocatable :: ind_var_act_species(:)       !< Indices of variable activity species in the chemical system (primary then secondary variable activity)
    contains
    !> Set procedures
        procedure :: set_ind_mins_chem_syst                  !< Set mineral indices in chemical system
        procedure :: set_ind_non_flow_species                 !< Set all non-flowing species
        procedure :: set_single_ind_non_flow_species          !< Set single non-flowing species by index
        procedure :: set_num_non_flow_species             !< Set number of non-flowing species
        procedure :: set_num_solids                          !< Set number of solid phases
        procedure :: set_chem_syst_react_zone                !< Associate chemical system pointer
        procedure :: set_cat_exch_zone                       !< Set cation exchange zone
        procedure :: set_ind_eq_reacts                       !< Set equilibrium reaction indices
        procedure :: set_stoich_mat_react_zone               !< Set stoichiometric matrix
        procedure :: set_ind_mins_stoich_mat                 !< Set mineral indices in stoich matrix
        procedure :: set_ind_gases_stoich_mat                !< Set gas indices in stoich matrix
        procedure :: set_num_mins_cst_act                    !< Set number of constant activity minerals
        procedure :: set_num_mins_var_act                    !< Set number of variable activity minerals
        procedure :: set_num_mins                            !< Set total number of minerals
        procedure :: set_gas_phase                           !< Set gas phase object
        procedure :: set_speciation_alg_dimensions           !< Set speciation algebra dimensions
        procedure :: set_CV_params                           !< Set convergence parameters
        procedure :: set_ind_var_act_species                 !< Set aqueous var act species indices
    !> Allocate/deallocate procedures
        procedure :: allocate_ind_non_flow_species            !< Allocate non-flowing species array
        procedure :: allocate_ind_mins                       !< Allocate mineral indices array
        procedure :: allocate_ind_eq_reacts                  !< Allocate equilibrium reaction indices
        procedure :: allocate_ind_gases_stoich_mat           !< Allocate gas indices in stoich matrix
        procedure :: allocate_ind_var_act_species            !< Allocate indices of variable activity species in chemical system
        procedure :: deallocate_react_zone                   !< Deallocate all allocatable components
    !> Update procedures
        procedure :: update_num_eq_reacts                    !< Update number of equilibrium reactions
        procedure :: update_eq_reactions                     !< Update equilibrium reactions
    !> Get procedures
        procedure :: get_eq_csts_react_zone                  !< Get equilibrium constants array of reactive zone
        procedure :: get_Se_nc_react_zone                     !< Get stoichiometric matrix of secondary variable activity species in reactive zone
        procedure :: get_Se2v_react_zone                     !< Get stoichiometric matrix of secondary variable activity species in reactive zone
    !> Compute procedures
        procedure :: compute_U_SkT_prod                      !< Compute U*S_k^T product
        procedure :: compute_speciation_alg_arrays           !< Compute speciation algebra arrays
    !> Query procedures
        procedure :: is_nf_species_in_react_zone             !< Check if non-flowing species in zone
        procedure :: is_mineral_in_react_zone                !< Check if mineral in zone
    !> Copy procedures
        procedure :: copy_react_zone                         !< Copy from another reactive zone
    !> Rearrange procedures
        procedure :: rearrange_ind_non_flow_species           !< Rearrange non-flowing species in standard order
    end type
!**************************************************************************************************
    !> \brief External procedure interfaces for reactive zone operations
    !> \details
    !>   These interfaces define external procedures that operate on reactive zones:
    !>   - Stoichiometric matrix setup (solid, gas, full)
    !>   - I/O operations (read, write)
    !>   - Reaction updates and comparisons
    interface
         !> \brief Set stoichiometric matrix for reactive zone
         !> \details
         !>   Constructs the full stoichiometric matrix for all equilibrium reactions
         !>   in the reactive zone. Matrix rows are species, columns are reactions.
         !> \param[in,out] this Reactive zone object
         subroutine set_stoich_mat_react_zone(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
         end subroutine
        
         !> \brief Set solid stoichiometric matrix
         !> \details
         !>   Constructs stoichiometric matrix for solid phase reactions only.
         !> \param[in,out] this Reactive zone object
         subroutine set_stoich_mat_sol_rz(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
         end subroutine
         
         !> \brief Set gas stoichiometric matrix
         !> \details
         !>   Constructs stoichiometric matrix for gas phase reactions only.
         !> \param[in,out] this Reactive zone object
         subroutine set_stoich_mat_gas_rz(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
        end subroutine
        
       
        
        
        

        
       
        
       
        
        !> \brief Read reactive zone data for Lagrangian approach
        !> \details
        !>   Reads reactive zone configuration from file for Lagrangian transport.
        !>   Includes non-flowing species, equilibrium reactions, and zone properties.
        !> \param[in,out] this Reactive zone object
        !> \param[in] filename Input file name
        !> \param[in] line Line number to start reading
        subroutine read_reactive_zone_Lagr(this,filename,line)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
            character(len=*), intent(in) :: filename !> nombre del archivo de entrada
            integer(kind=4), intent(in) :: line
        end subroutine
        
        !> \brief Write reactive zone to output
        !> \details
        !>   Writes complete reactive zone information including:
        !>   - Zone ID
        !>   - Non-flowing species
        !>   - Equilibrium reactions
        !>   - Stoichiometric matrices
        !> \param[in] this Reactive zone object
        subroutine write_reactive_zone(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
        end subroutine
        
       
        

        
       
        
      
        
        !> \brief Update equilibrium reactions after zone changes
        !> \details
        !>   Updates the equilibrium reactions in reactive zone based on
        !>   changes to non-flowing species or reaction indices.
        !> \param[in,out] this Reactive zone object
        !> \param[in] old_eq_reacts_ind Indices of old equilibrium reactions
        subroutine update_eq_reactions(this,old_eq_reacts_ind)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: old_eq_reacts_ind(:)
        end subroutine
        
        !> \brief Compare two reactive zones for equivalence
        !> \details
        !>   Checks if two reactive zones have the same non-flowing species.
        !>   Used to detect zone changes in reactive transport.
        !> \param[in] react_zone_1 First reactive zone
        !> \param[in] react_zone_2 Second reactive zone
        !> \param[out] flag True if same non-flowing species, false otherwise
        subroutine compare_react_zones(react_zone_1,react_zone_2,flag)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c), intent(in) :: react_zone_1
            class(reactive_zone_c), intent(in) :: react_zone_2
            logical, intent(out) :: flag !> true if same non flowing species, false otherwise
        end subroutine
    end interface
    
    
    
    contains
    !******************************************************************************************************************************
        !> \brief Set number of non-flowing species
        !> \details
        !>   Sets the count of immobile species in the reactive zone.
        !>   If not provided, computes as: num_minerals + num_surf_compl + num_gases_eq
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_non_flow_species Optional: explicit count
        subroutine set_num_non_flow_species(this,num_non_flow_species)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_non_flow_species
            if (present(num_non_flow_species)) then
                this%num_non_flow_species=num_non_flow_species
            else
                this%num_non_flow_species=this%num_minerals+this%cat_exch_zone%num_surf_compl+this%gas_phase%num_gases_eq
            end if
        end subroutine
    !******************************************************************************************************************************
        !> \brief Set number of solid phases
        !> \details
        !>   Sets total count of solid phases (minerals + surface complexes).
        !>   If not provided, computes from component counts.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_solids Optional: explicit count
        subroutine set_num_solids(this,num_solids)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_solids
            if (present(num_solids)) then
                this%num_solids=num_solids
            else
                this%num_solids=this%num_minerals+this%cat_exch_zone%num_surf_compl
            end if
        end subroutine
    !******************************************************************************************************************************
        !> \brief Set number of constant activity minerals
        !> \details
        !>   Sets count of minerals with constant activity (infinite reservoir).
        !>   These minerals buffer solution composition.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_mins_cst_act Number of constant activity minerals
        subroutine set_num_mins_cst_act(this,num_mins_cst_act)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_cst_act
            this%num_minerals_cst_act=num_mins_cst_act
        end subroutine
    !******************************************************************************************************************************
        !> \brief Set number of variable activity minerals
        !> \details
        !>   Sets count of minerals with variable activity (finite amount).
        !>   These minerals can dissolve/precipitate completely.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_mins_var_act Number of variable activity minerals
        subroutine set_num_mins_var_act(this,num_mins_var_act)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_var_act
            this%num_minerals_var_act=num_mins_var_act
        end subroutine
    !******************************************************************************************************************************
        !> \brief Set gas phase object
        !> \details
        !>   Assigns the gas phase component of the reactive zone.
        !>   Gas phase may contain constant or variable activity gases.
        !> \param[in,out] this Reactive zone object
        !> \param[in] gas_phase Gas phase object to assign
        subroutine set_gas_phase(this,gas_phase)
            implicit none
            class(reactive_zone_c) :: this
            class(gas_phase_c), intent(in) :: gas_phase
            this%gas_phase=gas_phase
        end subroutine
    !******************************************************************************************************************************
        !> \brief Check consistency of solid phase count
        !> \details
        !>   Verifies that num_solids equals num_minerals + num_surf_compl.
        !>   Stops execution if inconsistency detected.
        !>   
        !>   Validation equation: num_solids = num_minerals + num_surf_compl
        !> \param[in] this Reactive zone object (const)
        subroutine check_num_solids(this)
            implicit none 
            class(reactive_zone_c) :: this !< Reactive zone to validate
            !< Check if total solid count matches sum of components
            if (this%num_solids/=this%num_minerals+this%cat_exch_zone%num_surf_compl) then
                error stop "Wrong number of solids in reactive zone object" !< Terminate if inconsistent
            end if
        end subroutine
    !******************************************************************************************************************************
        !> \brief Allocate non-flowing species array
        !> \details
        !>   Allocates memory for non-flowing species array.
        !>   If num_non_flow_species not provided, computes from components.
        !>   Deallocates existing array if already allocated.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_non_flow_species Optional: explicit count
        subroutine allocate_ind_non_flow_species(this,num_non_flow_species)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_non_flow_species
            if (present(num_non_flow_species)) then
                call this%set_num_non_flow_species(num_non_flow_species)
            else
                call this%set_num_non_flow_species()
            end if
            if (allocated(this%ind_non_flow_species)) then
                deallocate(this%ind_non_flow_species)
            end if
            allocate(this%ind_non_flow_species(this%num_non_flow_species))
        end subroutine

        
    !******************************************************************************************************************************
        !> \brief Set indices non-flowing species array
        !> \details
        !>   Sets the complete non-flowing species indices array for the reactive zone.
        !>   
        !>   If ind_non_flow_species provided: direct assignment
        !>   If not provided: constructs from chemical system in same order:
        !>   1. Free surface
        !>   2. Variable activity minerals in equilibrium
        !>   3. Surface complexes
        !>   4. Variable activity gases in equilibrium
        !>   5. Constant activity minerals in equilibrium
        !>   6. Constant activity gases in equilibrium
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_non_flow_species Optional: explicit species array
        subroutine set_ind_non_flow_species(this,ind_non_flow_species)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_non_flow_species(:)
            
            integer(kind=4) :: i,counter
            
            if (present(ind_non_flow_species)) then
                this%ind_non_flow_species=ind_non_flow_species
            else
                ! j=0 !> counter constant non flowing species
                ! k=0 !> counter variable non flowing species
                counter=0
                call this%allocate_ind_non_flow_species(this%num_solids+this%gas_phase%num_gases_eq)
                !> First: surface complexes
                do i=1,this%cat_exch_zone%num_surf_compl
                    this%ind_non_flow_species(counter+1)=this%chem_syst%speciation_alg%num_aq_prim_species+1
                    this%ind_non_flow_species(counter+this%num_minerals_var_act+i)=this%chem_syst%speciation_alg%num_prim_species+&
                        this%chem_syst%speciation_alg%num_aq_sec_var_act_species+&
                        this%chem_syst%num_minerals_eq_var_act+i    !> Indexing in chem system
                end do
                counter=counter+this%cat_exch_zone%num_surf_compl-this%cat_exch_zone%num_exch_cats
                !> Second: minerals with variable activity
                do i=1,this%num_minerals_var_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_prim_species+&
                        this%chem_syst%speciation_alg%num_aq_sec_var_act_species+i    !> Indexing in chem system
                end do
                counter=counter+this%num_minerals_var_act+this%cat_exch_zone%num_exch_cats
                !> Third: gases with variable activity
                do i=1,this%gas_phase%num_gases_eq_var_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_prim_species+&
                        this%chem_syst%speciation_alg%num_aq_sec_var_act_species+&
                        this%chem_syst%num_minerals_eq_var_act+&
                        this%chem_syst%cat_exch_zone%num_exch_cats+i    !> Indexing in chem system
                end do
                counter=counter+this%num_minerals_var_act
                !> Fourth: minerals with constant activity
                do i=1,this%num_minerals_cst_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-&
                        this%chem_syst%gas_phase%num_cst_act_species-&
                        this%chem_syst%num_minerals_eq_cst_act+i    !> Indexing in chem system
                end do
                counter=counter+this%num_minerals_cst_act
                !> Fifth: gases with constant activity
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-&
                        this%chem_syst%gas_phase%num_gases_eq_cst_act+i    !> Indexing in chem system
                end do
                counter=counter+this%gas_phase%num_gases_eq_cst_act
                !!> Fifth: gases with variable activity
                !do i=1,this%gas_phase%num_gases_eq_var_act
                !    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-&
                !        this%chem_syst%gas_phase%num_gases_eq_cst_act-&
                !        this%chem_syst%num_minerals_eq_cst_act-&
                !        this%chem_syst%num_minerals_kin_cst_act-&
                !        this%chem_syst%aq_phase%wat_flag-&
                !        this%chem_syst%gas_phase%num_gases_eq_var_act+i    !> Indexing in chem system
                !end do
                !counter=counter+this%gas_phase%num_gases_eq_var_act
            end if
            !print *, this%ind_non_flow_species
        end subroutine
    !******************************************************************************************************************************
        !> \brief Rearrange non-flowing species in standard order
        !> \details
        !>   Reorganizes ind_non_flow_species array into canonical order:
        !>   1. Constant activity minerals
        !>   2. Constant activity gases
        !>   3. Variable activity minerals
        !>   4. Surface complexes
        !>   5. Variable activity gases
        !>   
        !>   This ordering is required for consistent indexing in reactive transport.
        !> \param[in,out] this Reactive zone object
        subroutine rearrange_ind_non_flow_species(this)
            implicit none
            class(reactive_zone_c) :: this
            
            integer(kind=4) :: i,j,k,l,m,n,counter
            !type(species_c), allocatable :: old_non_flow_species(:)
            LOGICAL :: flag_gas,flag_surf
            
            counter=0
                !> First: minerals with constant activity
                do i=1,this%num_minerals_cst_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-this%chem_syst%gas_phase%num_gases_eq_cst_act-&
                        this%chem_syst%num_minerals_eq_cst_act+i    !> Indexing in chem system
                end do
                counter=counter+this%num_minerals_cst_act
                !> Second: gases with constant activity
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-&
                        this%chem_syst%gas_phase%num_gases_eq_cst_act+i    !> Indexing in chem system
                end do
                counter=counter+this%gas_phase%num_gases_eq_cst_act
                !> Third: minerals with variable activity
                do i=1,this%num_minerals_var_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_prim_species+&
                        this%chem_syst%speciation_alg%num_aq_sec_var_act_species+&
                        this%chem_syst%num_minerals_kin_var_act+i    !> Indexing in chem system
                end do
                counter=counter+this%num_minerals_var_act
                !> Fourth: surface complexes
                do i=1,this%cat_exch_zone%num_surf_compl
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_prim_species+&
                        this%chem_syst%speciation_alg%num_aq_sec_var_act_species+&
                        this%chem_syst%num_minerals_kin_var_act+&
                        this%chem_syst%num_minerals_eq_var_act+i    !> Indexing in chem system
                end do
                counter=counter+this%cat_exch_zone%num_surf_compl
                !> Fifth: gases with variable activity
                do i=1,this%gas_phase%num_gases_eq_var_act
                    this%ind_non_flow_species(counter+i)=this%chem_syst%speciation_alg%num_species-&
                        this%chem_syst%gas_phase%num_gases_eq_cst_act-&
                        this%chem_syst%num_minerals_eq_cst_act-&
                        this%chem_syst%num_minerals_kin_cst_act-&
                        this%chem_syst%aq_phase%wat_flag-&
                        this%chem_syst%gas_phase%num_gases_eq_var_act+i    !> Indexing in chem system
                end do
                counter=counter+this%gas_phase%num_gases_eq_var_act
            !old_non_flow_species=this%ind_non_flow_species
            !deallocate(this%ind_non_flow_species)
            !call this%allocate_ind_non_flow_species()
            
            !j=0 !> counter constant activity gases
            !k=0 !> counter variable activity gases
            !l=0 !> counter surface complexes
            !m=0 !> counter constant activity minerals
            !n=0 !> counter variable activity minerals
            !do i=1,this%num_non_flow_species
            !    call old_non_flow_species(I)%is_gas(flag_gas)
            !    if (flag_gas .eqv. .true.) then
            !        if (old_non_flow_species(i)%cst_act_flag.eqv..true.) then
            !            j=j+1
            !            call this%ind_non_flow_species(this%num_minerals_cst_act+j)%copy_species(old_non_flow_species(I))
            !        else
            !            k=k+1
            !            call this%ind_non_flow_species(&
            !            this%num_non_flow_species-this%gas_phase%num_gases_eq_var_act+k)%copy_species(&
            !            old_non_flow_species(I))
            !        end if
            !    else
            !        call old_non_flow_species(I)%is_surf_compl(flag_surf)
            !        if (flag_surf .eqv. .true.) then
            !            l=l+1
            !            call this%ind_non_flow_species(this%num_minerals+this%gas_phase%num_gases_eq_cst_act+l)%copy_species(&
            !            old_non_flow_species(I))
            !        else if (old_non_flow_species(i)%cst_act_flag .eqv. .true.) then
            !            m=m+1
            !            call this%ind_non_flow_species(m)%copy_species(old_non_flow_species(I))
            !        else
            !            n=n+1
            !        call this%ind_non_flow_species(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act+n)%copy_species(&
            !            old_non_flow_species(I))
            !        end if
            !    end if
            !end do
        end subroutine
    !******************************************************************************************************************************
        !> \brief Set single non-flowing species by index
        !> \details
        !>   Sets one non-flowing species at specified position using species
        !>   from chemical system.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_non_flow_species_ind Index in reactive zone (1-based)
        !> \param[in] chem_syst_ind Index in chemical system species array
        subroutine set_single_ind_non_flow_species(this,ind_rz,ind_chem_syst)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: ind_rz !> index in reactive zone
            integer(kind=4), intent(in) :: ind_chem_syst !> index in chemical system
            if (ind_rz<1 .or. ind_rz>ind_chem_syst .or. ind_chem_syst>this%chem_syst%speciation_alg%num_species) then
                error stop "Index out of bounds"
            else
                this%ind_non_flow_species(ind_rz)=ind_chem_syst
            end if
        end subroutine
    !******************************************************************************************************************************
        !> \brief Associate chemical system with reactive zone
        !> \details
        !>   Sets pointer to chemical system object.
        !>   Chemical system must remain in scope for lifetime of reactive zone.
        !> \param[in,out] this Reactive zone object
        !> \param[in] chem_syst Chemical system to associate (target must persist)
        subroutine set_chem_syst_react_zone(this,chem_syst)
            implicit none
            class(reactive_zone_c) :: this
            class(chem_system_c), intent(in), target :: chem_syst
            this%chem_syst=>chem_syst
        end subroutine
    !******************************************************************************************************************************
        !> \brief Get equilibrium constants for reactive zone
        !> \details
        !>   Retrieves equilibrium constants for all reactions in proper order.
        !>   For gas reactions, can adjust K by constant activity gases:
        !>   K_adj = K_eq / P_gas (for constant pressure gases)
        !>   
        !>   Ordering depends on speciation_alg%flag_comp:
        !>   - If TRUE: cst act minerals, cst act gases, redox, aq complex, var act minerals, exch, var act gases
        !>   - If FALSE: redox, aq complex, all minerals, exch, all gases
        !> \param[in] this Reactive zone object (const)
        !> \param[in] cst_act_gases Optional: constant activity gas pressures
        !> \return K Array of equilibrium constants
        function get_eq_csts_react_zone(this,cst_act_gases) result(K)
            implicit none
            class(reactive_zone_c), intent(in) :: this !> reactive zone object
            real(kind=8), intent(in), optional :: cst_act_gases(:) !> constant activity gases in gas phase
            real(kind=8), allocatable :: K(:) !> equilibrium constants
            
            integer(kind=4) :: i,ind_eq_react
          
            allocate(K(this%speciation_alg%num_eq_reactions))
            
            if (this%speciation_alg%num_eq_reactions/=size(this%ind_eq_reacts)) error stop "Wrong size of ind_eq_reacts array"

            ind_eq_react=0 !> counter equilibrium reactions

            if (this%speciation_alg%flag_comp .eqv. .true.) then !> If component matrix has NO variable activity species
                do i=1,this%num_minerals_cst_act
                    K(i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(i))%eq_cst
                end do
                ind_eq_react=ind_eq_react+this%num_minerals_cst_act
                if (present(cst_act_gases)) then
                    do i=1,this%gas_phase%num_gases_eq_cst_act
                        K(this%num_minerals_cst_act+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                            this%num_minerals_cst_act+i))%eq_cst/cst_act_gases(i)
                    end do
                else
                    do i=1,this%gas_phase%num_gases_eq_cst_act
                        K(this%num_minerals_cst_act+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                            this%num_minerals_cst_act+i))%eq_cst
                    end do
                end if
                ind_eq_react=ind_eq_react+this%gas_phase%num_gases_eq_cst_act
                do i=ind_eq_react+1,this%speciation_alg%num_eq_reactions
                    K(i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(i))%eq_cst
                end do
            else !> If component matrix has variable activity species
                do i=1,this%chem_syst%num_aq_eq_reacts
                    K(i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(i))%eq_cst
                end do
                ind_eq_react=ind_eq_react+this%chem_syst%num_aq_eq_reacts !> counter equilibrium reactions
                do i=1,this%num_minerals
                    K(ind_eq_react+i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(ind_eq_react+i))%eq_cst
                end do
                ind_eq_react=ind_eq_react+this%num_minerals
                do i=1,this%cat_exch_zone%num_exch_cats
                    K(ind_eq_react+i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(ind_eq_react+i))%eq_cst
                end do
                ind_eq_react=ind_eq_react+this%cat_exch_zone%num_exch_cats
                if (present(cst_act_gases)) then
                    do i=1,this%gas_phase%num_gases_eq
                        K(ind_eq_react+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                            ind_eq_react+i))%eq_cst/cst_act_gases(i)
                    end do
                else
                    do i=1,this%gas_phase%num_gases_eq
                        K(ind_eq_react+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                            ind_eq_react+i))%eq_cst
                    end do
                end if
            end if
            
            
           
        end function
        
     
       
        
        
        
      
       
        
        !> \brief Set equilibrium reaction indices in reactive zone
        !> \details
        !>   Assigns indices of equilibrium reactions from the chemical system,
        !>   ordering them according to speciation algebra requirements.
        !>   
        !>   **Ordering depends on flag_comp:**
        !>   
        !>   If flag_comp = TRUE (component matrix has NO variable activity species):
        !>   1. Constant activity minerals
        !>   2. Constant activity gases
        !>   3. Redox reactions
        !>   4. Aqueous complexation
        !>   5. Variable activity minerals
        !>   6. Cation exchange
        !>   7. Variable activity gases
        !>   
        !>   If flag_comp = FALSE (component matrix HAS variable activity species):
        !>   1. Redox reactions
        !>   2. Aqueous complexation
        !>   3. Variable activity minerals
        !>   4. Constant activity minerals
        !>   5. Cation exchange
        !>   6. Variable activity gases
        !>   7. Constant activity gases
        !>   
        !>   This ordering is critical for proper matrix construction in speciation.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_eq_reacts Optional: explicit indices array (must match num_eq_reactions)
        subroutine set_ind_eq_reacts(this,ind_eq_reacts)
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_eq_reacts(:) !> indices of equilibrium reactions in chemical system
            
            integer(kind=4) :: i,l,j,k,eq_react_ind,n_eq,sp_ind  !< Loop indices and counters
            integer(kind=4), allocatable :: eq_react_indices(:)  !< Temporary indices array
            logical :: flag  !< Flag for species search results

            call this%allocate_ind_eq_reacts()  !< Allocate based on num_eq_reactions
            
            !> \subsection explicit_indices Case 1: Explicit indices provided
            if (present(ind_eq_reacts)) then
                !call this%allocate_eq_reactions(size(eq_reactions_ind))
                if (size(ind_eq_reacts)/=this%speciation_alg%num_eq_reactions) error stop "Error: Size mismatch in equilibrium reaction indices" !< Validate size matches expected
                this%ind_eq_reacts=ind_eq_reacts !< Direct assignment
                
            !> \subsection both_comp_true Case 2: Both reactive zone AND chemical system use component matrix WITH variable activity species
            else if ((this%speciation_alg%flag_comp .eqv. .true.) .and. (this%chem_syst%speciation_alg%flag_comp .eqv. .true.)) then
                !> Step 2a: Assign aqueous equilibrium reactions (redox + complexation)
                !> These go in middle: after cst_act minerals/gases, before var_act species
                do i=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act+i)=&
                        this%chem_syst%num_minerals_eq_cst_act+this%chem_syst%gas_phase%num_gases_eq_cst_act+i
                end do
                
                !> Step 2b: Initialize counters for non-flowing species assignment
                i=1 !> Counter: constant activity minerals & gases in reactive zone
                l=1 !> Counter: variable activity minerals, surface complexes & gases in reactive zone
                j=1 !> Counter: equilibrium reactions in chemical system (search index)
                k=1 !> Counter: non-flowing species in reactive zone
                
                !> Step 2c: Loop through non-flowing species to find their reactions
                if (this%num_non_flow_species>0) then
                    do !< Infinite loop with exit conditions
                        !> Skip placeholder species 'x-' (temporary fix)
                        if (this%chem_syst%species(this%ind_non_flow_species(k))%name=='x-') then
                            k=k+1 !< Move to next non-flowing species
                        end if
                        
                        !> Search for equilibrium reaction containing this non-flowing species
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(&
                            this%ind_non_flow_species(k),flag,sp_ind)
                        
                        if (flag .eqv. .true.) then !< Found reaction for this species
                            !> Assign index based on constant vs variable activity
                            if (this%chem_syst%species(this%ind_non_flow_species(k))%cst_act_flag .eqv. .true.) then !< Constant activity species
                                !> Place at beginning of ind_eq_reacts array
                                this%ind_eq_reacts(i)=j
                                i=i+1 !< Increment constant activity counter
                            else !< Variable activity species
                                !> Place after: cst_act mins + cst_act gases + aqueous reactions
                                this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act + & 
                                    this%chem_syst%num_aq_eq_reacts+l)=j
                                l=l+1 !< Increment variable activity counter
                            end if
                            
                            !> Move to next non-flowing species or exit if done
                            if (k<this%num_non_flow_species) then
                                k=k+1 !< Next non-flowing species
                                j=1   !< Reset reaction search to beginning
                            else
                                exit !< All non-flowing species processed
                            end if
                        else if (j<this%chem_syst%speciation_alg%num_eq_reactions) then !< Not found yet, continue searching
                            j=j+1 !< Try next equilibrium reaction
                        else !< Searched all reactions without finding one for this species
                            print *, this%chem_syst%species(this%ind_non_flow_species(k))%name
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
                
            !> \subsection rz_comp_true_cs_comp_false Case 3: Reactive zone uses component matrix WITH var act, chem system WITHOUT var act
            else if ((this%speciation_alg%flag_comp .eqv. .true.) .and. (this%chem_syst%speciation_alg%flag_comp .eqv. .false.))then
                !> Step 3a: Assign aqueous equilibrium reactions sequentially at beginning
                do i=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(i)=i !< Sequential 1:1 mapping
                end do
                
                !> Step 3b: Initialize counters for non-flowing species assignment
                i=1 !> Counter: constant activity minerals & gases in reactive zone
                l=1 !> Counter: variable activity minerals, surface complexes & gases in reactive zone
                j=1 !> Counter: equilibrium reactions in chemical system (search index)
                k=1 !> Counter: non-flowing species in reactive zone
                
                !> Step 3c: Loop through non-flowing species to find their reactions
                if (this%num_non_flow_species>0) then
                    do !< Infinite loop with exit conditions
                        !> Skip placeholder species 'x-' (temporary fix)
                        if (this%chem_syst%species(this%ind_non_flow_species(k))%name=='x-') then
                            k=k+1 !< Move to next non-flowing species
                        end if
                        
                        !> Search for equilibrium reaction containing this non-flowing species
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(&
                            this%ind_non_flow_species(k),flag,sp_ind)
                        
                        if (flag .eqv. .true.) then !< Found reaction for this species
                            !> Assign index based on constant vs variable activity
                            if (this%chem_syst%species(this%ind_non_flow_species(k))%cst_act_flag .eqv. .true.) then !< Constant activity species
                                !> Place after aqueous reactions
                                this%ind_eq_reacts(this%chem_syst%num_aq_eq_reacts+i)=j
                                i=i+1 !< Increment constant activity counter
                            else !< Variable activity species
                                !> Place after: cst_act mins + cst_act gases + aqueous reactions
                                this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act + & 
                                    this%chem_syst%num_aq_eq_reacts+l)=j
                                l=l+1 !< Increment variable activity counter
                            end if
                            
                            !> Move to next non-flowing species or exit if done
                            if (k<this%num_non_flow_species) then
                                k=k+1 !< Next non-flowing species
                                j=1   !< Reset reaction search to beginning
                            else
                                exit !< All non-flowing species processed
                            end if
                        else if (j<this%chem_syst%speciation_alg%num_eq_reactions) then !< Not found yet, continue searching
                            j=j+1 !< Try next equilibrium reaction
                        else !< Searched all reactions without finding one for this species
                            print *, this%chem_syst%species(this%ind_non_flow_species(k))%name
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
                
            !> \subsection default_case Case 4: Default case (flag_comp = FALSE for reactive zone)
            !> Simpler ordering: aqueous reactions first, then all non-flowing species reactions
            else
                !> Step 4a: Initialize counters
                i=0 !> Counter: equilibrium reactions in reactive zone
                j=1 !> Counter: equilibrium reactions in chemical system (search index)
                k=1 !> Counter: non-flowing species in reactive zone
                
                !> Step 4b: Assign all aqueous equilibrium reactions sequentially
                do l=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(l)=l !< Sequential 1:1 mapping
                end do
                i=i+this%chem_syst%num_aq_eq_reacts+1 !< Move index past aqueous reactions
                
                !> Step 4c: Loop through non-flowing species to find their reactions
                if (this%num_non_flow_species>0) then
                    do !< Infinite loop with exit conditions
                        !> Skip placeholder species 'x-' (temporary fix)
                        if (this%chem_syst%species(this%ind_non_flow_species(k))%name=='x-') then
                            k=k+1 !< Move to next non-flowing species
                        end if
                        
                        !> Search for equilibrium reaction containing this non-flowing species
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(this%ind_non_flow_species(k),flag,sp_ind)
                        
                        if (flag .eqv. .true.) then !< Found reaction for this species
                            this%ind_eq_reacts(i)=j !< Assign index sequentially (no cst/var distinction)
                            i=i+1 !< Increment reaction counter
                            
                            !> Move to next non-flowing species or exit if done
                            if (k<this%num_non_flow_species) then
                                k=k+1 !< Next non-flowing species
                                j=1   !< Reset reaction search to beginning
                            else
                                exit !< All non-flowing species processed
                            end if
                        else if (j<this%chem_syst%speciation_alg%num_eq_reactions) then !< Not found yet, continue searching
                            j=j+1 !< Try next equilibrium reaction
                        else !< Searched all reactions without finding one for this species
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
            end if
        end subroutine
                
        !> \brief Set mineral indices in chemical system
        !> \details
        !>   Maps minerals in reactive zone to their indices in the chemical system.
        !>   
        !>   If ind_mins_chem_syst provided: direct assignment
        !>   If not provided: defaults to all equilibrium minerals in chem system
        !>   
        !>   Indices refer to positions in chem_syst%minerals array.
        !>   Validates that indices don't exceed available minerals.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_mins_chem_syst Optional: mineral indices in chemical system
        subroutine set_ind_mins_chem_syst(this,ind_mins_chem_syst)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_mins_chem_syst(:) !> indices of minerals in chemical system
            integer(kind=4) :: i  !< Loop index
            
            if (present(ind_mins_chem_syst)) then  !< Explicit indices provided
                if (size(ind_mins_chem_syst)>this%chem_syst%num_minerals_eq) error stop  !< Validate count
                this%num_minerals=size(ind_mins_chem_syst)  !< Set count from array size
                this%ind_mins_chem_syst=ind_mins_chem_syst  !< Direct assignment
            else !> default: all minerals in chemical system
                this%num_minerals=this%chem_syst%num_minerals_eq  !< Use all equilibrium minerals
                do i=1,this%num_minerals  !< Assign indices sequentially
                    !< Skip kinetic minerals (first num_minerals_kin), take equilibrium minerals
                    this%ind_mins_chem_syst(i)=this%chem_syst%num_minerals_kin+i
                end do
                ! do i=1,this%num_minerals
                !     if (this%minerals(i)%mineral%cst_act_flag.eqv..true.) then
                !         this%num_minerals_cst_act=this%num_minerals_cst_act+1
                !     else
                !         this%num_minerals_var_act=this%num_minerals_var_act+1
                !     end if
                ! end do
            end if
        end subroutine
        
        !> \brief Allocate mineral indices arrays
        !> \details
        !>   Allocates memory for both ind_mins_chem_syst and ind_mins_stoich_mat.
        !>   Deallocates existing arrays if already allocated.
        !>   
        !>   Both arrays have size num_minerals.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_minerals Optional: explicit mineral count (validated against chem system)
        subroutine allocate_ind_mins(this,num_minerals)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_minerals !> number of minerals in reactive zone
            
            if (present(num_minerals)) then  !< Set mineral count if provided
                if (num_minerals<0 .or. num_minerals>this%chem_syst%num_minerals_eq) error stop  !< Validate range
                this%num_minerals=num_minerals
            end if
            !< Deallocate if already allocated
            if (allocated(this%ind_mins_chem_syst)) deallocate(this%ind_mins_chem_syst)
            if (allocated(this%ind_mins_stoich_mat)) deallocate(this%ind_mins_stoich_mat)
            !< Allocate both indices arrays with same size
            allocate(this%ind_mins_chem_syst(this%num_minerals),this%ind_mins_stoich_mat(this%num_minerals))
        end subroutine
        
       
        
        !> \brief Update number of equilibrium reactions
        !> \details
        !>   Decrements num_eq_reactions by number of old reactions.
        !>   Used when updating reactive zone after removing equilibrium reactions.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_old_eq_reacts Number of reactions to remove
        subroutine update_num_eq_reacts(this,num_old_eq_reacts)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_old_eq_reacts  !< Count of reactions being removed
            this%speciation_alg%num_eq_reactions=this%speciation_alg%num_eq_reactions-num_old_eq_reacts  !< Decrement total
        end subroutine
        
        !> \brief Allocate equilibrium reaction indices array
        !> \details
        !>   Allocates ind_eq_reacts array based on num_eq_reactions in speciation_alg.
        !>   Deallocates existing array if already allocated.
        !> \param[in,out] this Reactive zone object
        subroutine allocate_ind_eq_reacts(this)
            implicit none
            class(reactive_zone_c) :: this
            if (allocated(this%ind_eq_reacts)) then  !< Deallocate if exists
                deallocate(this%ind_eq_reacts)
            end if
            !< Allocate with size from speciation algebra
            allocate(this%ind_eq_reacts(this%speciation_alg%num_eq_reactions))
        end subroutine
        
      
        
        !> \brief Check if non-flowing species is in reactive zone
        !> \details
        !>   Searches for a non-flowing species by name in the reactive zone.
        !>   Returns flag indicating presence and optionally the index.
        !>   
        !>   Used for checking mineral, gas, or surface complex presence.
        !> \param[in] this Reactive zone object (const)
        !> \param[in] nf_species Non-flowing species to search for
        !> \param[out] flag TRUE if species found, FALSE otherwise
        !> \param[out] nf_species_ind Index in ind_non_flow_species array (0 if not found)
        subroutine is_nf_species_in_react_zone(this,ind_sp_chem_syst,flag,ind_sp_rz)
            implicit none
            class(reactive_zone_c), intent(in) :: this
            integer(kind=4), intent(in) :: ind_sp_chem_syst  !< index in chemical system species array
            logical, intent(out) :: flag  !< TRUE if found
            integer(kind=4), intent(out) :: ind_sp_rz  !< Index if found, 0 otherwise
            
            integer(kind=4) :: i  !< Loop index

            ind_sp_rz=0  !< Initialize to "not found"
            flag=.false.  !< Initialize to "not found"
            
            do i=1,this%num_non_flow_species  !< Loop through all non-flowing species
                if (ind_sp_chem_syst==this%ind_non_flow_species(i)) then  !< Compare names
                    flag=.true.  !< Mark as found
                    ind_sp_rz=i  !< Store index
                    exit  !< Stop searching
                end if
            end do
        end subroutine
        
       
        
       
        
        
        
        !> \brief Check if mineral is in reactive zone
        !> \details
        !>   Searches for a mineral in the reactive zone by name and activity flag.
        !>   Both name AND constant activity flag must match.
        !>   
        !>   Searches through minerals referenced in ind_mins_chem_syst.
        !> \param[in] this Reactive zone object (const)
        !> \param[in] mineral Mineral to search for
        !> \param[out] flag TRUE if mineral found with same activity flag, FALSE otherwise
        !> \param[out] index Optional: index in reactive zone minerals (0 if not found)
        subroutine is_mineral_in_react_zone(this,ind_min_chem_syst,flag,ind_min_rz)
            class(reactive_zone_c), intent(in) :: this !> reactive zone
            integer(kind=4), intent(in) :: ind_min_chem_syst !> index in chemical system minerals array
            logical, intent(out) :: flag !> true if mineral is in reactive zone, false otherwise
            integer(kind=4), intent(out) :: ind_min_rz !> index of mineral in reactive zone
            
            integer(kind=4) :: i  !< Loop index
            
            flag=.false.  !< Initialize to "not found"
            ind_min_rz=0  !< Initialize to "not found"
            !if (present(ind_min_rz)) then
            !    index=0  !< Initialize index to "not found"
            !end if
            
            do i=1,this%num_minerals  !< Loop through all minerals in reactive zone
                !< Check both name and constant activity flag match
                if (ind_min_chem_syst==this%ind_mins_chem_syst(i)) then
                    flag=.true.  !< Mark as found
                    ind_min_rz=i  !< Store index
                    exit  !< Stop searching
                end if
            end do
        end subroutine
        
       
        
        !> \brief Set cation exchange zone
        !> \details
        !>   Assigns the cation exchange component of the reactive zone.
        !>   
        !>   If cat_exch_zone provided: validates all surface complexes exist in
        !>   chemical system, then assigns.
        !>   
        !>   If not provided: defaults to chemical system's cation exchange zone.
        !>   
        !>   Stops execution if any surface complex not found in chemical system.
        !> \param[in,out] this Reactive zone object
        !> \param[in] cat_exch_zone Optional: explicit cation exchange zone
        subroutine set_cat_exch_zone(this,cat_exch_zone)
            implicit none
            class(reactive_zone_c) :: this
            class(cat_exch_zone_c), intent(in), optional :: cat_exch_zone
            
            integer(kind=4) :: i  !< Loop index for surface complexes
            logical :: flag  !< Flag for surface complex validation
            
            if (present(cat_exch_zone)) then  !< Explicit cation exchange zone provided
                !< Validate all surface complexes exist in chemical system
                do i=1,cat_exch_zone%num_surf_compl
                    call this%chem_syst%cat_exch_zone%is_surf_compl_in(cat_exch_zone%surf_compl(i),flag)
                    if (flag.eqv..false.) then  !< Surface complex not in chem system
                        error stop "Surface complex not in chemical system"
                    end if
                end do
                this%cat_exch_zone=cat_exch_zone  !< Assign after validation
            else  !< No cation exchange zone provided
                this%cat_exch_zone=this%chem_syst%cat_exch_zone  !< Use chemical system's default
            end if
        end subroutine
        
        !> \brief Deallocate reactive zone components
        !> \details
        !>   Frees all allocatable arrays in the reactive zone.
        !>   Must be called before object goes out of scope to prevent memory leaks.
        !>   
        !>   Deallocates:
        !>   - ind_mins_chem_syst
        !>   - ind_non_flow_species
        !>   - stoich_mat
        !>   - ind_mins_stoich_mat
        !>   - ind_gases_stoich_mat
        !>   - ind_eq_reacts
        !>   - U_SkT_prod
        !> \param[in,out] this Reactive zone object
        subroutine deallocate_react_zone(this)
            implicit none
            class(reactive_zone_c) :: this
            deallocate(this%ind_mins_chem_syst)  !< Free mineral indices in chem system
            deallocate(this%ind_non_flow_species)  !< Free non-flowing species array
            deallocate(this%stoich_mat)  !< Free stoichiometric matrix
            deallocate(this%ind_mins_stoich_mat)  !< Free mineral indices in stoich matrix
            deallocate(this%ind_gases_stoich_mat)  !< Free gas indices in stoich matrix
            deallocate(this%ind_eq_reacts)  !< Free equilibrium reaction indices
            deallocate(this%U_SkT_prod)  !< Free U*S_k^T product matrix
        end subroutine
        
        !> \brief Assign reactive zone from another instance (chapuza)
        !> \details
        !>   Copies all components from source reactive zone to this one.
        !>   Performs deep copy of allocatable arrays and pointer associations.
        !>   
        !>   **Critical:** Pointers (chem_syst, CV_params) must be associated in source.
        !>   
        !>   Copies:
        !>   - Pointers (chem_syst, CV_params)
        !>   - All allocatable arrays
        !>   - All scalar values
        !>   - Embedded objects (gas_phase, cat_exch_zone, speciation_alg)
        !>   
        !>   Stops execution if required pointers not associated.
        !> \param[in,out] this Reactive zone object (target)
        !> \param[in] react_zone Reactive zone to copy from (source)
        subroutine copy_react_zone(this,react_zone) !> chapuza 
            implicit none
            class(reactive_zone_c) :: this !< reactive zone
            class(reactive_zone_c), intent(in) :: react_zone !< reactive zone to be assigned
            
            !< Check and associate chemical system pointer
            if (associated(react_zone%chem_syst)) then
                this%chem_syst=>react_zone%chem_syst  !< Associate pointer (shallow copy)
            else
                error stop "Chemical system not associated with reactive zone"
            end if
            
            !< Check and associate convergence parameters pointer
            if (associated(react_zone%CV_params)) then
                this%CV_params=>react_zone%CV_params  !< Associate pointer (shallow copy)
            else
                error stop "CV parameters not associated with reactive zone"
            end if
            
            !< Copy mineral indices if allocated
            if (allocated(react_zone%ind_mins_chem_syst)) then
                this%ind_mins_chem_syst=react_zone%ind_mins_chem_syst  !< Deep copy array
            end if
            
            !< Copy non-flowing species if allocated
            if (allocated(react_zone%ind_non_flow_species)) then
                this%ind_non_flow_species=react_zone%ind_non_flow_species  !< Deep copy array
            end if
            
            this%num_non_flow_species=react_zone%num_non_flow_species  !< Copy count
            this%gas_phase=react_zone%gas_phase  !< Copy gas phase object
            this%num_minerals=react_zone%num_minerals  !< Copy mineral count
            this%num_minerals_cst_Act=react_zone%num_minerals_cst_Act  !< Copy cst act mineral count
            this%num_minerals_var_Act=react_zone%num_minerals_var_Act  !< Copy var act mineral count
            this%cat_exch_zone=react_zone%cat_exch_zone  !< Copy cation exchange zone
            this%num_solids=react_zone%num_solids  !< Copy solid count
            
            !< Copy stoichiometric matrix if allocated
            if (allocated(react_zone%stoich_mat)) then
                this%stoich_mat=react_zone%stoich_mat  !< Deep copy matrix
            end if
            
            !< Copy mineral indices in stoich matrix if allocated
            if (allocated(react_zone%ind_mins_stoich_mat)) then
                this%ind_mins_stoich_mat=react_zone%ind_mins_stoich_mat  !< Deep copy array
            end if
            
            !< Copy gas indices in stoich matrix if allocated
            if (allocated(react_zone%ind_gases_stoich_mat)) then
                this%ind_gases_stoich_mat=react_zone%ind_gases_stoich_mat  !< Deep copy array
            end if
            
            this%speciation_alg=react_zone%speciation_alg  !< Copy speciation algebra object
            
            !< Copy variable activity species indices if allocated
            if (allocated(react_zone%ind_var_act_species)) then
                if (allocated(this%ind_var_act_species)) deallocate(this%ind_var_act_species)
                this%ind_var_act_species=react_zone%ind_var_act_species  !< Deep copy array
            end if
            
            !< Copy equilibrium reaction indices if allocated
            if (allocated(react_zone%ind_eq_reacts)) then
                call this%set_ind_eq_reacts(react_zone%ind_eq_reacts)  !< Use setter for proper assignment
            end if
            
            !< Copy U*S_k^T product if kinetic reactions exist
            if (react_zone%chem_syst%num_kin_reacts>0 .and. allocated(react_zone%U_SkT_prod)) then
                this%U_SkT_prod=react_zone%U_SkT_prod  !< Deep copy matrix
            end if
        end subroutine
        
        !> \brief Set speciation algebra object
        !> \details
        !>   Assigns the speciation algebra component of the reactive zone.
        !>   Deep copy of speciation_algebra_s object.
        !> \param[in,out] this Reactive zone object
        !> \param[in] speciation_alg Speciation algebra object to assign
        subroutine set_speciation_alg(this,speciation_alg)
            implicit none
            class(reactive_zone_c) :: this !< Reactive zone to modify
            type(speciation_algebra_s), intent(in) :: speciation_alg !< Speciation algebra to assign
            this%speciation_alg=speciation_alg !< Deep copy assignment
        end subroutine
        
        !> \brief Set speciation algebra dimensions based on reactive zone composition
        !> \details
        !>   Computes and sets the dimensions for speciation algebra object.
        !>   
        !>   **Calculated dimensions:**
        !>   - n_sp: Total species (aqueous + non-flowing)
        !>   - n_c: Constant activity species count
        !>   - n_eq: Total equilibrium reactions
        !>   - flag_cat_exch: TRUE if cation exchange present
        !>   
        !>   **Constant activity species:**
        !>   - Water (if present in aqueous phase)
        !>   - Constant activity minerals
        !>   - Constant activity gases
        !>   
        !>   **Equilibrium reactions:**
        !>   - Aqueous complexation
        !>   - Redox reactions
        !>   - Mineral dissolution/precipitation
        !>   - Gas dissolution
        !>   - Cation exchange
        !>   
        !>   If no non-flowing species: uses only aqueous chemistry from chem_syst.
        !> \param[in,out] this Reactive zone object
        !> \param[in] flag_comp Optional: TRUE if component matrix has no constant activity species
        subroutine set_speciation_alg_dimensions(this,flag_comp)
            implicit none
            class(reactive_zone_c) :: this
            logical, intent(in), optional :: flag_comp !> TRUE if component matrix has no constant activity species (De Simoni et al, 2005), FALSE otherwise
            
            integer(kind=4) :: i,n_sp,n_c,n_eq,n_gas_kin  !< Dimension counters
            logical :: flag_cat_exch  !< Flag for cation exchange presence
            
            n_gas_kin=0  !< Initialize kinetic gas counter to zero
            
            !> \subsection validate_chem_syst Validate chemical system association
            if (.not. associated(this%chem_syst)) then  !< Chemical system pointer must be set
                error stop "Chemical system not associated with reactive zone"
                
            !> \subsection case_nf_species Case 1: Reactive zone has non-flowing species (minerals, gases, surf complexes)
            else if (this%num_non_flow_species>0) then
                !> Compute total number of species in reactive zone
                !> n_sp = aqueous species + non-flowing species
                !> Note: kinetic minerals excluded (they don't participate in equilibrium)
                n_sp=this%chem_syst%aq_phase%num_species+this%num_non_flow_species!+this%chem_syst%num_minerals_kin
                
                !> Compute number of constant activity species
                !> Includes:
                !> - Water (if present, via wat_flag which is 0 or 1)
                !> - Constant activity minerals (pure phases buffering solution)
                !> - Constant activity gases (partial pressure fixed)
                n_c=this%chem_syst%aq_phase%wat_flag+this%num_minerals_cst_act+this%gas_phase%num_cst_act_species
                !do i=1,this%num_minerals  !< Alternative: loop through minerals checking cst_act_flag
                !    if (this%minerals(I)%mineral%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                !do i=1,this%chem_syst%num_minerals_kin  !< Alternative: check kinetic minerals
                !    if (this%chem_syst%minerals(i)%mineral%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                !do i=1,this%gas_phase%num_species  !< Alternative: loop through gases
                !    if (this%gas_phase%gases(i)%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                
                !> Compute total number of equilibrium reactions
                !> Includes:
                !> - Aqueous complexation reactions (aq species formation)
                !> - Redox reactions (in aq_eq_reacts)
                !> - Mineral dissolution/precipitation (num_minerals reactions)
                !> - Gas dissolution (num_gases_eq reactions)
                !> - Cation exchange (num_exch_cats reactions)
                n_eq=this%num_minerals+this%gas_phase%num_gases_eq+this%cat_exch_zone%num_exch_cats+ & 
                    this%chem_syst%num_aq_eq_reacts
                
                !> Determine if cation exchange is active
                if (this%cat_exch_zone%num_surf_compl>0) then
                    flag_cat_exch=.true.  !< Cation exchange present (surface complexation)
                else
                    flag_cat_exch=.false.  !< No cation exchange
                end if
                
            !> \subsection case_no_nf_species Case 2: No non-flowing species (purely aqueous chemistry)
            else !> all equilibrium reactions are aqueous (no non-flowing species)
                !> Use only aqueous phase dimensions from chemical system
                n_sp=this%chem_syst%aq_phase%num_species  !< Only aqueous species count
                n_c=this%chem_syst%aq_phase%num_cst_act_species  !< Only aqueous constant activity species (water, etc.)
                n_eq=this%chem_syst%num_aq_eq_reacts  !< Only aqueous equilibrium reactions
                
                !> Check for cation exchange in chemical system
                !> (even though no non-flowing species in reactive zone, chem system may have it)
                if (this%chem_syst%cat_exch_zone%num_surf_compl>0) then
                    flag_cat_exch=.true.  !< Cation exchange defined in chem system
                else
                    flag_cat_exch=.false.  !< No cation exchange
                end if
            end if
            
            !> \subsection set_flags Set optional flags
            !> Set component matrix flag if provided
            !> flag_comp = TRUE: component matrix U excludes constant activity species (De Simoni approach)
            !> flag_comp = FALSE: component matrix U includes constant activity species
            if (present(flag_comp)) then
                call this%speciation_alg%set_flag_comp(flag_comp)
            end if
            
            !> Set cation exchange flag in speciation algebra
            call this%speciation_alg%set_flag_cat_exch(flag_cat_exch)
            
            !> \subsection set_dimensions Set all speciation algebra dimensions
            !> Arguments:
            !> 1. n_sp: total species count
            !> 2. n_eq: total equilibrium reactions count
            !> 3. n_c: constant activity species count
            !> 4. num_aq_species: aqueous species count (for indexing)
            !> 5. num_aq_species - wat_flag: variable activity aqueous species count
            call this%speciation_alg%set_dimensions(n_sp,n_eq,n_c,this%chem_syst%aq_phase%num_species, & 
                this%chem_syst%aq_phase%num_species-this%chem_syst%aq_phase%wat_flag)
        end subroutine
        
        !> \brief Compute speciation algebra arrays and handle matrix modifications
        !> \details
        !>   Computes all speciation algebra arrays (component matrix, inverse Se, logK_tilde)
        !>   and handles necessary species/reaction swaps to maintain matrix invertibility.
        !>   
        !>   **Workflow:**
        !>   1. Get equilibrium constants K
        !>   2. Call speciation_alg%compute_arrays to compute U, inv(Se), logK_tilde
        !>   3. If swap needed (flag=TRUE): update indices accordingly
        !>   
        !>   **Swap handling:**
        !>   - If flag_comp=FALSE: swap reactions in ind_eq_reacts
        !>   - If flag_comp=TRUE: swap species in ind_var_act_species
        !> \param[in,out] this Reactive zone object
        !> \param[out] flag TRUE if stoichiometric matrix modified, FALSE otherwise
        !> \param[out] swap Species/reaction indices to swap [2] (already allocated)
        !> \param[in] cst_act_gases Optional: constant activity gas pressures
        subroutine compute_speciation_alg_arrays(this,flag,swap,cst_act_gases)
            implicit none
            class(reactive_zone_c) :: this !< Reactive zone object
            logical, intent(out) :: flag !< TRUE if stoichiometric matrix has been modified, FALSE otherwise
            integer(kind=4), intent(out) :: swap(:) !< Species or reactions to swap in stoichiometric matrix (already allocated)
            real(kind=8), intent(in), optional :: cst_act_gases(:) !< Constant activity of gases in reactive zone (chapuza)

            !< Local variable declarations
            real(kind=8), allocatable :: Se(:,:) !< Stoichiometric matrix (temporary)
            real(kind=8), allocatable :: K(:) !< Equilibrium constants array
            real(kind=8), allocatable :: aux_Se(:,:) !< Auxiliary stoichiometric matrix for swapping
            real(kind=8), allocatable :: aux_Sk(:,:) !< Auxiliary kinetic stoichiometric matrix
            integer(kind=4) :: aux_col !< Auxiliary column index
            type(species_c), allocatable :: aux_species(:) !< Auxiliary species array for swapping
            !type(eq_reaction_c), allocatable :: aux_eq_reacts(:) !< Auxiliary equilibrium reactions
            
            !logical :: flag
            !type(aq_phase_c), target :: aux_aq_phase
                        
            !call aq_phase_new%copy_attributes(this%aq_phase)

            flag=.false. !< Initialize to FALSE: no modification by default
            swap=0 !< Initialize to zero: no swap by default
            
            !> \par Case 1: Reactive zone has equilibrium reactions
            if (this%speciation_alg%num_eq_reactions>0) then
                !Se=this%stoich_mat
                K=this%get_eq_csts_react_zone(cst_act_gases) !< Get equilibrium constants
                !< Compute component matrix, inverse Se, logK_tilde
                call this%speciation_alg%compute_arrays(this%stoich_mat,K,this%CV_params%zero,flag,swap)
                !< Check if matrix modification (swap) was required
                if (flag .eqv. .true.) then
                    !allocate(aux_Se(this%speciation_Alg%num_eq_reactions,this%speciation_alg%num_Species))
                    !aux_Se=Se
                    !> \par Swap Case A: Component matrix WITHOUT variable activity species
                    if (this%speciation_alg%flag_comp .eqv. .false.) then
                        !allocate(aux_ind_eq_reacts(2))
                        !aux_ind_eq_reacts=this%ind_eq_reacts(swap)
                        !< Swap reaction indices to maintain matrix properties
                        this%ind_eq_reacts(swap(1))=swap(2) !< Assign second reaction to first position
                        this%ind_eq_reacts(swap(2))=swap(1) !< Assign first reaction to second position
                        !this%stoich_mat(swap(1),:)=aux_Se(swap(2),:)
                        !this%stoich_mat(swap(2),:)=aux_Se(swap(1),:)
                        !deallocate(aux_ind_eq_reacts)
                    !> \par Swap Case B: Component matrix WITH variable activity species
                    else
                        !< Allocate auxiliary arrays for species swap
                        allocate(aux_Sk(this%chem_syst%speciation_alg%num_species,this%chem_syst%speciation_alg%num_species),aux_species(2))
                        aux_species=this%chem_syst%species(swap) !< Store species to swap
                        !this%chem_syst%species(swap(1))=aux_species(2)
                        !this%chem_syst%species(swap(2))=aux_species(1)
                        aux_Sk=this%chem_syst%Sk !< Store kinetic stoichiometric matrix
                        !this%chem_syst%Sk(:,swap(1))=aux_Sk(:,swap(2))
                        !this%chem_syst%Sk(:,swap(2))=aux_Sk(:,swap(1))
                        ! this%stoich_mat(:,swap(1))=aux_Se(:,swap(2))
                        ! this%stoich_mat(:,swap(2))=aux_Se(:,swap(1))
                        deallocate(aux_species,aux_Sk) !< Free temporary arrays
                        !< Swap variable activity species indices
                        this%ind_var_act_species(swap(1))=swap(2) !< Assign second species to first position
                        this%ind_var_act_species(swap(2))=swap(1) !< Assign first species to second position
                        !print *, "Swapped species ", swap(1), " and ", swap(2) !< Debug output
                    end if
                    !deallocate(aux_Se)
                end if
                deallocate(K) !< Free equilibrium constants array
            !> \par Case 2: No equilibrium reactions in reactive zone, use chemical system
            else if (associated(this%chem_syst)) then
                Se=this%chem_syst%Se !< Get stoichiometric matrix from chemical system
                K=this%chem_syst%get_eq_csts() !< Get equilibrium constants from chemical system
                !> \par Case 2A: Component matrix includes constant activity species
                if (this%speciation_alg%flag_comp .eqv. .false.) then
                    !< Set component matrix for constant activity species as identity
                    this%speciation_alg%comp_mat_cst_act=id_matrix(this%speciation_alg%num_prim_species)
                    !< Compute inverse of secondary species stoichiometric matrix
                    call this%speciation_alg%compute_inv_Se_2(Se(&
                        :,this%speciation_alg%num_prim_species+1:this%speciation_alg%num_species),&
                        this%CV_params%zero)
                    !< Compute modified equilibrium constants
                    call this%speciation_alg%compute_logK_tilde(K)
                !> \par Case 2B: Component matrix excludes constant activity species
                else
                    !< Set component matrix as identity (variable activity only)
                    this%speciation_alg%comp_mat=id_matrix(this%speciation_alg%num_prim_species)
                    !< Set aqueous component matrix as identity
                    this%speciation_alg%comp_mat_aq=id_matrix(this%speciation_alg%num_aq_prim_species)
                    !< Compute inverse of non-constant activity secondary species stoichiometric matrix
                    !call this%speciation_alg%compute_inv_Se_nc_2(Se(&
                    !    :,this%speciation_alg%num_prim_species+1:this%speciation_alg%num_var_act_species),&
                    !    this%CV_params%zero)
                    !call this%speciation_alg%compute_logK_tilde(K)
                end if
                !call this%speciation_alg%compute_comp_mat()
                !call this%speciation_alg%compute_arrays(this%chem_syst%Se,K,this%CV_params%zero,flag,swap)
            !> \par Case 3: No chemical system associated - error condition
            else
                error stop "Chemical system not associated with reactive zone" !< Fatal error
            end if
        end subroutine
        
        !> \brief Set convergence parameters pointer
        !> \details
        !>   Associates the CV_params pointer with provided convergence parameters.
        !>   Target must remain in scope for lifetime of reactive zone.
        !> \param[in,out] this Reactive zone object
        !> \param[in] CV_params Convergence parameters object (target must persist)
        subroutine set_CV_params(this,CV_params)
            implicit none
            class(reactive_zone_c) :: this
            class(CV_params_s), intent(in), target :: CV_params  !< Convergence parameters (must be target)
            this%CV_params=>CV_params  !< Associate pointer (shallow copy)
        end subroutine
        
        !> \brief Compute U*S_k^T product matrix for kinetic reactions
        !> \details
        !>   Computes the product of component matrix U and transpose of kinetic
        !>   stoichiometric matrix S_k. This product is used in reactive transport
        !>   to transform kinetic reaction rates to primary component space.
        !>   
        !>   **Matrix dimensions:**
        !>   - If num_eq_reactions == 0: U_SkT_prod = S_k^T (identity transformation)
        !>     Size: [num_var_act_species × num_kin_reacts]
        !>   - If num_eq_reactions > 0: U_SkT_prod = U * S_k^T
        !>     Size: [num_prim_species × num_kin_reacts]
        !>   
        !>   **Kinetic reaction selection:**
        !>   - If ind_kin provided: uses only specified kinetic reactions
        !>   - If not provided: uses all kinetic reactions in chemical system
        !>   
        !>   The product transforms: dc/dt = U * S_k^T * r_k
        !>   where r_k are kinetic reaction rates.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_kin Optional: indices of kinetic reactions in chemical system
        subroutine compute_U_SkT_prod(this,ind_kin)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_kin(:) !> indices of kinetic reactions in chemical system
            
            !> \subsection case_specific_kin Case A: Specific kinetic reactions provided via ind_kin
            if (present(ind_kin)) then  !< Use only specified kinetic reactions
                
                !> Case A1: No equilibrium reactions → Identity transformation (U = I)
                if (this%speciation_alg%num_eq_reactions==0) then
                    !< Allocate matrix if needed
                    !< Dimensions: [num_var_act_species × number of specified kinetic reactions]
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_var_act_species,size(ind_kin)))
                    end if
                    !< Direct assignment: U_SkT_prod = S_k^T (transpose of kinetic stoichiometric matrix)
                    !< Extract rows for specified kinetic reactions (ind_kin)
                    !< Extract columns for variable activity species (1:num_var_act_species)
                    !< Transpose to get [species × reactions]
                    this%U_SkT_prod=transpose(this%chem_syst%Sk(ind_kin,1:this%speciation_alg%num_var_act_species)) !> chapuza
                    
                !> Case A2: Equilibrium reactions present → Component matrix transformation (U * S_k^T)
                else
                    !< Allocate matrix if needed
                    !< Dimensions: [num_prim_species × total num_kin_reacts in chem_syst]
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,size(ind_kin)))
                    end if
                    !< Matrix multiplication: U_SkT_prod = U * S_k^T
                    !< U: component matrix [num_prim × num_var_act]
                    !< S_k^T: transposed kinetic stoich [num_var_act × num_kin_specified]
                    !< Result: [num_prim × num_kin_specified]
                    this%U_SkT_prod = matmul(this%speciation_alg%comp_mat, &
                                    transpose(this%chem_syst%Sk(ind_kin, 1:this%speciation_alg%num_var_act_species)))
                end if
                
            !> \subsection case_all_kin Case B: Use all kinetic reactions in chemical system (default)
            else !> we consider all kinetic reactions in the chemical system
                
                !> Case B1: No equilibrium reactions → Identity transformation
                if (this%speciation_alg%num_eq_reactions==0) then
                    !< Allocate matrix if needed
                    !< Dimensions: [num_var_act_species × total num_kin_reacts in chem_syst]
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_var_act_species,this%chem_syst%num_kin_reacts))
                    end if
                    !< Direct assignment: U_SkT_prod = S_k^T (all kinetic reactions)
                    !< Extract all rows of S_k (all kinetic reactions)
                    !< Extract columns for variable activity species only
                    !< Transpose to get [species × reactions]
                    this%U_SkT_prod=transpose(this%chem_syst%Sk(:,1:this%speciation_alg%num_var_act_species)) !> chapuza
                    
                !> Case B2: Both equilibrium and kinetic reactions present → Component transformation
                else if (this%chem_syst%num_kin_reacts>0) then
                    !< Allocate matrix if needed
                    !< Dimensions: [num_prim_species × total num_kin_reacts in chem_syst]
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,this%chem_syst%num_kin_reacts))
                    end if
                    !< Matrix multiplication: U_SkT_prod = U * S_k^T (all kinetic reactions)
                    !< U: component matrix [num_prim × num_var_act]
                    !< S_k^T: transposed full kinetic stoich [num_var_act × num_kin_reacts]
                    !< Result: [num_prim × num_kin_reacts]
                    !< This transforms kinetic rates from concentration to component basis
                    this%U_SkT_prod = matmul(this%speciation_alg%comp_mat, &
                                             transpose(this%chem_syst%Sk(:, 1:this%speciation_alg%num_var_act_species)))
                    
                !> Case B3: No kinetic reactions in system → Empty matrix
                else if (.not. allocated(this%U_SkT_prod)) then
                    !< Allocate zero-width matrix (no kinetic reactions to transform)
                    !< Dimensions: [num_prim_species × 0]
                    allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,0))
                    this%U_SkT_prod=0d0  !< Initialize to zero (though unused)
                end if
            end if
        end subroutine
        
        !> \brief Set indices of variable activity species
        !> \details
        !>   Sets the indices of variable activity species in the chemical system.
        !>   
        !>   **Default ordering:**
        !>   1. Primary species (first)
        !>   2. Secondary variable activity species (second)
        !>   
        !>   **Assumption:** Species in chemical system are already ordered by
        !>   variable vs constant activity.
        !>   
        !>   By default, assigns sequential indices 1, 2, 3, ... n_var_act_species.
        !> \param[in,out] this Reactive zone object
        subroutine set_ind_var_act_species(this,ind_var_act_species) !> sets the indices of variable activity species in the chemical system
        !! first primary species, then secondary variable activity species
        !! we assume that the species in the chemical system are already ordered in variable & constant activity
            class(reactive_zone_c) :: this !> reactive zone
            integer(kind=4), intent(in), optional :: ind_var_act_species(:) !> indices of variable activity species
            integer(kind=4) :: i  !< Loop index

            if (present(ind_var_act_species)) then
                !< Assign provided indices
                this%ind_var_act_species = ind_var_act_species
            else
                !< Assign sequential indices (default ordering)
                do i=1,this%speciation_alg%num_var_act_species
                    this%ind_var_act_species(i)=i !> by default: index equals position
                end do
            end if
        end subroutine

        !> \brief Set mineral indices in stoichiometric matrix
        !> \details
        !>   Maps minerals to their column positions in the stoichiometric matrix.
        !>   
        !>   **Default ordering if not provided:**
        !>   1. Variable activity minerals (after prim + sec aq species)
        !>   2. Constant activity minerals (after var act mins + exch cats + var act gases)
        !>   
        !>   **Matrix organization:**
        !>   - Columns: [primary aq | secondary aq | var act mins | exch cats | var act gases | water | cst act mins | cst act gases]
        !>   
        !>   Validates that provided indices don't exceed num_solids.
        !>   Assumes ind_mins_chem_syst is already allocated.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_mins_stoich_mat Optional: explicit indices in stoich matrix
        subroutine set_ind_mins_stoich_mat(this,ind_mins_stoich_mat)
        !> This subroutine sets the "ind_mins_stoich_mat" attribute
            !! We assume that the indices of minerals are already allocated
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_mins_stoich_mat(:) !> indices of solids in stoichiometric matrix
            integer(kind=4) :: i,num_mins,num_sp  !< Loop index, mineral counter, species counter
            
            if (present(ind_mins_stoich_mat)) then  !< Explicit indices provided
                if (size(ind_mins_stoich_mat)>this%num_solids) error stop  !< Validate size
                !this%num_solids=size(ind_mins_stoich_mat)
                this%ind_mins_stoich_mat=ind_mins_stoich_mat  !< Direct assignment
            else  !< Compute default indices
                num_mins=0 !> counter minerals
                num_sp=this%speciation_alg%num_prim_species+this%speciation_alg%num_sec_aq_species !> counter species in stoichiometric matrix
                
                !> Variable activity minerals in equilibrium (first)
                do i=1,this%num_minerals_var_act
                    this%ind_mins_stoich_mat(num_mins+i)=num_sp+i  !< Sequential indices after aq species
                end do
                num_mins=num_mins+this%num_minerals_var_act  !< Update mineral counter
                !< Skip over var act mins, exch cats, var act gases
                num_sp=num_sp+this%num_minerals_var_act+this%cat_exch_zone%num_exch_cats+this%gas_phase%num_gases_eq_var_act
                ! !> Surface complexes
                ! do i=1,this%cat_exch_zone%num_exch_cats
                !     this%ind_solids(num_mins+i)=num_sp+i
                ! end do
                ! num_mins=num_mins+this%cat_exch_zone%num_exch_cats
                ! num_sp=num_sp+this%cat_exch_zone%num_exch_cats+this%gas_phase%num_gases_eq_var_act+this%chem_syst%aq_phase%wat_flag
                
                !> Constant activity minerals in equilibrium (second)
                do i=1,this%num_minerals_cst_act
                    this%ind_mins_stoich_mat(num_mins+i)=num_sp+i  !< Indices after var act species
                end do
            end if
        end subroutine

        !> \brief Set gas indices in stoichiometric matrix
        !> \details
        !>   Maps gases to their column positions in the stoichiometric matrix.
        !>   
        !>   **Default ordering if not provided:**
        !>   1. Variable activity gases (after prim + sec aq + var act mins + exch cats)
        !>   2. Constant activity gases (after var act gases + water + cst act mins)
        !>   
        !>   **Matrix organization:**
        !>   - Variable act gases: columns after minerals and exchange
        !>   - Constant act gases: columns near end of matrix
        !>   
        !>   Validates that provided indices don't exceed num_gases_eq.
        !>   Allocates ind_gases_stoich_mat array if needed.
        !> \param[in,out] this Reactive zone object
        !> \param[in] ind_gases_stoich_mat Optional: explicit indices in stoich matrix
        subroutine set_ind_gases_stoich_mat(this,ind_gases_stoich_mat)
        !> This subroutine sets the "ind_gases_stoich_mat" attribute
            !! We assume that the number of gases in the gas phase attribute is already set
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_gases_stoich_mat(:) !> indices of gases in stoichiometric matrix
            integer(kind=4) :: i,num_gas,num_sp  !< Loop index, gas counter, species counter
            
            call this%allocate_ind_gases_stoich_mat()  !< Allocate array
            
            if (present(ind_gases_stoich_mat)) then  !< Explicit indices provided
                if (size(ind_gases_stoich_mat)>this%gas_phase%num_gases_eq) error stop  !< Validate size
                !this%num_gases=size(ind_gases)
                this%ind_gases_stoich_mat=ind_gases_stoich_mat  !< Direct assignment
            else  !< Compute default indices
                !this%num_gases=this%gas_phase%num_gases_eq
                num_gas=0 !> counter gases
                !< Start after primary, secondary aq, and var act mins
                num_sp=this%speciation_alg%num_var_act_species-this%gas_phase%num_gases_eq_var_act !> counter species in stoichiometric matrix
                
                !> Variable activity gases in equilibrium (first)
                do i=1,this%gas_phase%num_gases_eq_var_act
                    this%ind_gases_stoich_mat(num_gas+i)=num_sp+i  !< Sequential indices
                end do
                num_gas=num_gas+this%gas_phase%num_gases_eq_var_act  !< Update gas counter
                !< Skip to constant activity species section
                num_sp=num_sp+this%speciation_alg%num_cst_act_species-this%gas_phase%num_gases_eq_cst_act
                
                !> Constant activity gases in equilibrium (second)
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    this%ind_gases_stoich_mat(num_gas+i)=num_sp+i  !< Indices in cst act section
                end do
            end if
        end subroutine

        ! subroutine allocate_ind_mins_stoich_mat(this)
        ! !> This subroutine allocates the "ind_mins_stoich_mat" attribute
        !     class(reactive_zone_c) :: this
        !     if (allocated(this%ind_mins_stoich_mat)) then
        !         deallocate(this%ind_mins_stoich_mat)
        !     end if
        !     allocate(this%ind_mins_stoich_mat(this%num_solids))
        ! end subroutine

        !> \brief Allocate gas indices in stoichiometric matrix array
        !> \details
        !>   Allocates ind_gases_stoich_mat array based on num_gases_eq in gas_phase.
        !>   Deallocates existing array if already allocated.
        !> \param[in,out] this Reactive zone object
        subroutine allocate_ind_gases_stoich_mat(this)
        !> This subroutine allocates the "ind_gases_stoich_mat" attribute
            class(reactive_zone_c) :: this
            if (allocated(this%ind_gases_stoich_mat)) then  !< Deallocate if exists
                deallocate(this%ind_gases_stoich_mat)
            end if
            !< Allocate with size from gas phase
            allocate(this%ind_gases_stoich_mat(this%gas_phase%num_gases_eq))
        end subroutine

        !> \brief Set total number of minerals
        !> \details
        !>   Sets the count of all minerals (constant + variable activity) in reactive zone.
        !>   Validates that count is non-negative.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_mins Total number of minerals
        !> \brief Set total number of minerals
        !> \details
        !>   Sets the count of all minerals (constant + variable activity) in reactive zone.
        !>   Validates that count is non-negative.
        !> \param[in,out] this Reactive zone object
        !> \param[in] num_mins Total number of minerals
        subroutine set_num_mins(this,num_mins)
        class(reactive_zone_c) :: this !< Reactive zone to modify
        integer(kind=4), intent(in) :: num_mins !< Number of minerals in reactive zone
        this%num_minerals=num_mins  !< Assign mineral count
        end subroutine
        
        !> \brief Allocate variable activity species indices array
        !> \details
        !>   Allocates ind_var_act_species array based on num_var_act_species.
        !>   Deallocates existing array if already allocated.
        !> \param[in,out] this Reactive zone object
        subroutine allocate_ind_var_act_species(this)
        class(reactive_zone_c) :: this !< Reactive zone to modify
        !< Deallocate if already exists
        if (allocated(this%ind_var_act_species)) then
            deallocate(this%ind_var_act_species) !< Free existing array
        end if
        !< Allocate with size from speciation algebra
        allocate(this%ind_var_act_species(this%speciation_alg%num_var_act_species))
        end subroutine
        
        !> \brief Get non-constant activity stoichiometric matrix for reactive zone
        !> \details
        !>   Extracts the stoichiometric matrix columns for variable activity species.
        !>   
        !>   If ind_var_act_species provided: extracts only specified columns
        !>   If not provided: extracts all variable activity species columns
        !>   
        !>   Result matrix dimensions: [num_eq_reactions × num_var_act_species]
        !> \param[in] this Reactive zone object (const)
        !> \param[in] ind_var_act_species Optional: indices of variable activity species
        !> \return Se_nc Stoichiometric matrix for non-constant activity species
        function get_Se2v_react_zone(this,ind_sec_var_act_species) result(Se2v)
            class(reactive_zone_c) :: this !< Reactive zone object
            integer(kind=4), intent(in), optional :: ind_sec_var_act_species(:) !< Indices of variable activity species in chemical system
            real(kind=8), allocatable :: Se2v(:,:) !< Stoichiometric matrix subset

            integer(kind=4) :: i !< Loop index
        
            !< Allocate result matrix
            allocate(Se2v(this%speciation_alg%num_eq_reactions,this%speciation_alg%num_eq_reactions))
            !< Extract columns based on provided indices
            if (present(ind_sec_var_act_species)) then
                !< Loop through specified indices
                do i=1,size(ind_sec_var_act_species)
                    !< Extract column for each specified species
                    Se2v(:,i)=this%stoich_mat(:,this%ind_var_act_species(ind_sec_var_act_species(i)))
                end do
            else
                !< Extract all secondary variable activity species columns (sequential)
                Se2v=this%stoich_mat(:,this%speciation_alg%num_prim_species+1:this%speciation_alg%num_var_act_species)
            end if
        end function

        function get_Se_nc_react_zone(this,ind_var_act_species) result(Se_nc)
            class(reactive_zone_c) :: this !< Reactive zone object
            integer(kind=4), intent(in), optional :: ind_var_act_species(:) !< Indices of variable activity species in chemical system
            real(kind=8), allocatable :: Se_nc(:,:) !< Stoichiometric matrix subset

            integer(kind=4) :: i !< Loop index
        
            !< Allocate result matrix
            allocate(Se_nc(this%speciation_alg%num_eq_reactions,this%speciation_alg%num_var_act_species))
            !< Extract columns based on provided indices
            if (present(ind_var_act_species)) then
                !< Loop through specified indices
                do i=1,size(ind_var_act_species)
                    !< Extract column for each specified species
                    Se_nc(:,i)=this%stoich_mat(:,this%ind_var_act_species(ind_var_act_species(i)))
                end do
            else
                !< Extract all variable activity species columns (sequential)
                Se_nc=this%stoich_mat(:,1:this%speciation_alg%num_var_act_species)
            end if
        end function
end module reactive_zone_m