!> \file mineral_zone_m.f90
!> \brief Mineral zone module for spatial distribution of mineral assemblages
!> \details
!>   Defines mineral zones - spatial regions with specific mineral assemblages.
!>   Each zone contains a subset of minerals from the chemical system.
!>   Minerals can be classified as:
!>   - Equilibrium vs kinetic reactions
!>   - Constant vs variable activity
!>   
!>   Used in reactive transport to track spatial heterogeneity of mineral distribution.
!>   Critical for efficient computation - only active minerals tracked per zone.
!>
!> \author Jordi Petchamé-Guerrero
!> \date October 2025

module mineral_zone_m
    use mineral_m, only: mineral_c                                            !< Mineral phase class
    use chem_system_m, only: chem_system_c                                    !< Chemical system class
    implicit none
    save
    private !< Make all entities private by default

!**********************************************************************************************************************************!
    !> \brief Mineral zone class - spatial region with specific mineral assemblage
    !> \details
    !>   Represents a spatial zone characterized by its mineral composition.
    !>   
    !>   Key features:
    !>   - Links to global chemical system
    !>   - Stores indices of active minerals (not full mineral objects)
    !>   - Separates minerals by equilibrium/kinetic and constant/variable activity
    !>   - Enables efficient computation by tracking only relevant minerals
    !>   
    !>   Mineral ordering in ind_min_chem_syst:
    !>   1. Kinetic minerals with variable activity
    !>   2. Kinetic minerals with constant activity
    !>   3. Equilibrium minerals with variable activity
    !>   4. Equilibrium minerals with constant activity
    !>   
    !>   Applications:
    !>   - Spatial heterogeneity in reactive transport
    !>   - Tracking mineral precipitation/dissolution zones
    !>   - Adaptive reactive zone management
    type, public :: mineral_zone_c
        integer(kind=4) :: id=0                                               !< Unique identifier for this mineral zone
        integer(kind=4) :: num_minerals=0                                     !< Total number of minerals in this zone
        integer(kind=4) :: num_minerals_eq=0                                  !< Number of minerals assumed at equilibrium
        integer(kind=4) :: num_minerals_kin=0                                 !< Number of minerals with kinetic reactions
        integer(kind=4) :: num_minerals_kin_cst_act=0                         !< Kinetic minerals with constant activity
        integer(kind=4) :: num_minerals_kin_var_act=0                         !< Kinetic minerals with variable activity
        integer(kind=4), allocatable :: ind_min_chem_syst(:)                  !< Indices of minerals in chemical system (ordered: kin var act, kin cst act, eq var act, eq cst act)
        integer(kind=4), allocatable :: ind_min_Sk(:)                         !< Indices of minerals in kinetic stoichiometric matrix (ordered: kin var act, kin cst act)
        integer(kind=4) :: num_minerals_eq_cst_act=0                          !< Equilibrium minerals with constant activity
        integer(kind=4) :: num_minerals_eq_var_act=0                          !< Equilibrium minerals with variable activity
        class(chem_system_c), pointer :: chem_syst => null()                  !< Pointer to associated chemical system
    contains
        procedure :: set_num_mins_min_zone                            !< Set total number of minerals
        procedure :: set_num_mins_eq_min_zone                         !< Set number of equilibrium minerals
        procedure :: set_num_mins_kin_min_zone                        !< Set number of kinetic minerals
        procedure :: set_num_mins_kin_cst_act_min_zone                !< Set kinetic constant activity count
        procedure :: set_num_mins_eq_cst_act_min_zone                 !< Set equilibrium constant activity count
        procedure :: set_num_mins_eq_var_act_min_zone                 !< Set equilibrium variable activity count
        procedure :: set_num_mins_kin_var_act_min_zone                !< Set kinetic variable activity count
        procedure :: set_chem_syst_min_zone                           !< Associate chemical system
        procedure :: set_ind_min_Sk                                   !< Set indices in kinetic stoichiometric matrix
        procedure :: set_id                                           !< Set zone identifier
        !procedure :: allocate_minerals_eq                            !< [COMMENTED] Allocate equilibrium minerals array
        !procedure :: allocate_minerals_min_zone                      !< [COMMENTED] Allocate minerals array
        procedure :: allocate_ind_chem_syst_min_zone                  !< Allocate indices array for chemical system
        procedure :: allocate_ind_min_Sk                              !< Allocate indices array for kinetic matrix
        !procedure :: set_minerals_eq                                 !< [COMMENTED] Set equilibrium minerals
        !procedure :: update_mineral_zone                             !< [COMMENTED] Update zone configuration
        !procedure :: compare_min_zones                               !< [COMMENTED] Compare two mineral zones
        procedure :: is_mineral_in_min_zone                           !< Check if mineral belongs to zone
    end type
!**************************************************************************************************
    !> \brief Interface for updating mineral zone configuration
    !> \details
    !>   External interface for mineral zone update procedure.
    !>   Allows modification of mineral zone based on old mineral indices.
    !>   Used when mineral assemblage changes during simulation.
    interface
        !> \brief Update mineral zone based on old mineral indices
        !> \param[inout] this Mineral zone to update
        !> \param[in] old_min_ind Array of previous mineral indices
        subroutine update_mineral_zone(this,old_min_ind)
            import mineral_zone_c
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: old_min_ind(:)                     !< Previous mineral indices
        end subroutine
    end interface
    
    
    
    contains
    
        !> \brief Set total number of minerals in zone
        !> \details
        !>   Sets the total mineral count for this zone.
        !>   If num_minerals not provided, uses all minerals from chemical system.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_minerals Total number of minerals (optional)
        subroutine set_num_mins_min_zone(this,num_minerals)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in), optional :: num_minerals             !< Number of minerals
            if (present(num_minerals)) then                                   !< Check if optional argument was provided
                this%num_minerals=num_minerals                                !< Use provided value
            else                                                              !< If not provided
                this%num_minerals=this%chem_syst%num_minerals                 !< Use all minerals from chemical system
            end if                                                            !< End conditional block
        end subroutine
        
        !> \brief Set number of kinetic minerals
        !> \details
        !>   Sets count of minerals undergoing kinetic reactions.
        !>   If not provided, computed as: total - equilibrium minerals.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_min_kin Number of kinetic minerals (optional)
        subroutine set_num_mins_kin_min_zone(this,num_min_kin)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in), optional :: num_min_kin              !< Kinetic mineral count
            if (present(num_min_kin)) then                                    !< Check if optional argument was provided
                this%num_minerals_kin=num_min_kin                             !< Use provided value
            else                                                              !< If not provided
                this%num_minerals_kin=this%num_minerals-this%num_minerals_eq  !< Compute from total and equilibrium
            end if                                                            !< End conditional block
        end subroutine
        
        !> \brief Set number of equilibrium minerals
        !> \details
        !>   Sets count of minerals assumed to be at thermodynamic equilibrium.
        !>   Equilibrium minerals react instantaneously.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_min_eq Number of equilibrium minerals
        subroutine set_num_mins_eq_min_zone(this,num_min_eq)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: num_min_eq                         !< Equilibrium mineral count
            this%num_minerals_eq=num_min_eq                                  !< Assign equilibrium mineral count to member variable
        end subroutine
        
        !> \brief Set number of kinetic minerals with constant activity
        !> \details
        !>   Sets count of kinetic minerals with activity fixed at constant value.
        !>   Constant activity minerals don't affect solution speciation.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_mins_kin_cst_act Count of kinetic constant activity minerals
        subroutine set_num_mins_kin_cst_act_min_zone(this,num_mins_kin_cst_act)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: num_mins_kin_cst_act               !< Kinetic constant activity count
            this%num_minerals_kin_cst_act=num_mins_kin_cst_act               !< Assign kinetic constant activity count to member variable
        end subroutine
        
        !> \brief Set number of kinetic minerals with variable activity
        !> \details
        !>   Sets count of kinetic minerals with activity varying with composition.
        !>   Variable activity minerals affect solution speciation.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_mins_kin_var_act Count of kinetic variable activity minerals
        subroutine set_num_mins_kin_var_act_min_zone(this,num_mins_kin_var_act)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: num_mins_kin_var_act               !< Kinetic variable activity count
            this%num_minerals_kin_var_act=num_mins_kin_var_act               !< Assign kinetic variable activity count to member variable
        end subroutine

        !> \brief Set number of equilibrium minerals with constant activity
        !> \details
        !>   Sets count of equilibrium minerals with fixed activity.
        !>   Used for minerals at saturation with infinite reservoir.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_mins_eq_cst_act Count of equilibrium constant activity minerals
        subroutine set_num_mins_eq_cst_act_min_zone(this,num_mins_eq_cst_act)
            implicit none                                                     !< Disable implicit variable typing (implicit none missing in original)
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: num_mins_eq_cst_act                !< Equilibrium constant activity count
            this%num_minerals_eq_cst_act=num_mins_eq_cst_act                 !< Assign equilibrium constant activity count to member variable
        end subroutine
        
        !> \brief Associate chemical system with mineral zone
        !> \details
        !>   Sets pointer to the chemical system containing mineral definitions.
        !>   Zone references minerals from this system via indices.
        !> \param[inout] this Mineral zone object
        !> \param[in] chem_syst Chemical system to associate
        subroutine set_chem_syst_min_zone(this,chem_syst)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            class(chem_system_c), intent(in), target :: chem_syst             !< Chemical system (must be target)
            this%chem_syst=>chem_syst                                         !< Set pointer association to chemical system
        end subroutine
        
        !> \brief Allocate array for mineral indices in chemical system
        !> \details
        !>   Allocates ind_min_chem_syst array to store mineral indices.
        !>   Optionally sets num_minerals before allocation.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_mins Number of minerals (optional, sets num_minerals if provided)
        subroutine allocate_ind_chem_syst_min_zone(this,num_mins)
           implicit none                                                      !< Disable implicit variable typing
           class(mineral_zone_c) :: this                                      !< Mineral zone object
           integer(kind=4), intent(in), optional :: num_mins                  !< Number of minerals
           if (present(num_mins)) then                                        !< Check if optional argument was provided
               this%num_minerals=num_mins                                     !< Update mineral count
            end if                                                            !< End conditional block
            allocate(this%ind_min_chem_syst(this%num_minerals))               !< Allocate index array with size num_minerals
        end subroutine
        
        !> \brief Allocate array for kinetic mineral indices in stoichiometric matrix
        !> \details
        !>   Allocates ind_min_Sk array for indices in kinetic stoichiometric matrix.
        !>   Optionally sets num_minerals_kin before allocation.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_min_kin Number of kinetic minerals (optional)
        subroutine allocate_ind_min_Sk(this,num_min_kin)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in), optional :: num_min_kin              !< Number of kinetic minerals
            if (present(num_min_kin)) then                                    !< Check if optional argument was provided
                call this%set_num_mins_kin_min_zone(num_min_kin)              !< Update kinetic mineral count using setter method
            end if                                                            !< End conditional block
            !allocate(this%ind_min_kin_chem_syst(this%num_minerals_kin))     !< [COMMENTED] Alternative allocation (deprecated)
            allocate(this%ind_min_Sk(this%num_minerals_kin))                  !< Allocate kinetic index array with size num_minerals_kin
        end subroutine
        
      
        
        !> \brief Set equilibrium minerals (COMMENTED OUT)
        !> \details
        !>   Would set the equilibrium minerals for this zone.
        !>   Includes dimension validation.
        !>   Currently commented out - may be activated in future.
        !>
        !> \param[inout] this Mineral zone object
        !> \param[in] minerals_eq Array of equilibrium minerals
        !>
        !> Validates that size(minerals_eq) <= num_minerals
        
        !subroutine set_minerals_eq(this,minerals_eq)
        !    implicit none
        !    class(mineral_zone_c) :: this
        !    class(mineral_c), intent(in) :: minerals_eq(:)
        !    integer(kind=4) :: i,j
        !    
        !    if (size(minerals_eq)>this%num_minerals) error stop "Dimension error in minerals in equilibrium"
        !    this%minerals_eq=minerals_eq
        !end subroutine
        
        !> \brief Compare two mineral zones for equality (COMMENTED OUT)
        !> \details
        !>   Would compare two mineral zones to determine if they are identical.
        !>   
        !>   Comparison criteria:
        !>   - Same total number of minerals
        !>   - Same number of equilibrium minerals
        !>   - Same number of constant activity equilibrium minerals
        !>   - All mineral indices match between zones
        !>   
        !>   Currently commented out - may be used for zone optimization.
        !>
        !> \param[in] min_zone_1 First mineral zone
        !> \param[in] min_zone_2 Second mineral zone
        !> \param[out] flag TRUE if zones identical, FALSE otherwise
        !>
        !> Algorithm:
        !> 1. Check if mineral counts match
        !> 2. Check if equilibrium counts match
        !> 3. Check if constant activity counts match
        !> 4. Verify all mineral indices match
        
        ! subroutine compare_min_zones(min_zone_1,min_zone_2,flag)
        !     implicit none
        !     class(mineral_zone_c), intent(in) :: min_zone_1
        !     class(mineral_zone_c), intent(in) :: min_zone_2
        !     logical, intent(out) :: flag !> TRUE if same minerals & same properties, FALSE otherwise
            
        !     integer(kind=4) :: i,j
            
        !     flag=.true.
        !     if (min_zone_1%num_minerals/=min_zone_2%num_minerals) then
        !         flag=.false.
        !     else if (min_zone_1%num_minerals_eq/=min_zone_2%num_minerals_eq) then
        !         flag=.false.
        !     else if (min_zone_1%num_minerals_eq_cst_act /= min_zone_2%num_minerals_eq_cst_act) then
        !         flag = .false.
        !     else
        !         j=1 !> index of mineral in min_zone_2
        !         do i=1,min_zone_1%num_minerals
        !             if (min_zone_1%ind_chem_syst(i)/=min_zone_2%ind_chem_syst(j)) then
        !                 j=j+1
        !                 if (j>min_zone_2%num_minerals) then
        !                     flag=.false.
        !                     exit
        !                 end if
        !             else
        !                 j=1
        !             end if
        !         end do
        !     end if
        ! end subroutine
        
        !> \brief Allocate kinetic mineral index array (COMMENTED OUT)
        !> \details
        !>   Would allocate array for kinetic mineral indices.
        !>   Currently commented out - superseded by allocate_ind_min_Sk.
        !>
        !> \param[inout] this Mineral zone object
        
        ! subroutine allocate_ind_min_kin(this)
        !     implicit none
        !     class(mineral_zone_c) :: this
        !     allocate(this%ind_min_kin(this%num_minerals_kin))
        ! end subroutine
        
        !> \brief Check if mineral belongs to this zone
        !> \details
        !>   Searches for mineral in zone by comparing names.
        !>   Optionally returns the index in ind_min_chem_syst array.
        !>   
        !>   Search algorithm:
        !>   - Linear search through all minerals in zone
        !>   - Compares mineral names with chemical system minerals
        !>   - Returns on first match (early exit)
        !>
        !> \param[in] this Mineral zone to search
        !> \param[in] mineral Mineral to look for
        !> \param[out] flag TRUE if mineral found in zone, FALSE otherwise
        !> \param[out] mineral_ind Index in ind_min_chem_syst (0 if not found, optional)
        subroutine is_mineral_in_min_zone(this,mineral,flag,mineral_ind)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c), intent(in) :: this                         !< Mineral zone to search
            class(mineral_c), intent(in) :: mineral                           !< Mineral to find
            logical, intent(out) :: flag                                      !< TRUE if found
            integer(kind=4), intent(out), optional :: mineral_ind             !< Index in zone (0 if not found)
            
            integer(kind=4) :: i                                              !< Loop counter
            
            flag=.false.                                                      !< Initialize flag to "not found"
            if (present(mineral_ind)) then                                    !< Check if optional output argument was provided
                mineral_ind=0                                                 !< Initialize index to 0 (indicating not found)
            end if                                                            !< End conditional block
            do i=1,this%num_minerals                                          !< Loop through all minerals in zone
                !> Compare names with mineral from chemical system
                if (mineral%name==this%chem_syst%minerals(this%ind_min_chem_syst(i))%name) then  !< Check if names match
                    flag=.true.                                               !< Set flag to indicate mineral was found
                    if (present(mineral_ind)) then                            !< If output index was requested
                        mineral_ind=i                                         !< Return index in zone's mineral array
                    end if                                                    !< End conditional block
                    exit                                                      !< Early exit - mineral found, no need to continue searching
                end if                                                        !< End name comparison conditional
            end do                                                            !< End loop
        end subroutine

        !> \brief Set indices of kinetic minerals in stoichiometric matrix
        !> \details
        !>   Sets indices mapping kinetic minerals to kinetic stoichiometric matrix.
        !>   
        !>   Two modes of operation:
        !>   1. Explicit: Provide ind_min_Sk array directly
        !>   2. Automatic: Compute indices based on chemical system structure
        !>   
        !>   Automatic indexing algorithm:
        !>   - Variable activity kinetic minerals placed first
        !>   - Constant activity kinetic minerals placed after
        !>   - Accounts for primary species, surface complexes, gases, secondary species
        !>   
        !>   Index structure in kinetic stoichiometric matrix:
        !>   [primary species | kin var act mins | surf compl | gases | sec species | kin cst act mins]
        !>
        !> \param[inout] this Mineral zone object
        !> \param[in] ind_min_Sk Array of indices in kinetic matrix (optional)
        subroutine set_ind_min_Sk(this,ind_min_Sk)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in), optional :: ind_min_Sk(:)            !< Kinetic matrix indices
            integer(kind=4) :: i,num_min_kin,num_sp                           !< Loop counter, mineral counter, species offset
            
            if (present(ind_min_Sk)) then                                     !< Check if indices array was provided
                this%ind_min_Sk=ind_min_Sk                                    !< Use provided indices directly
            else                                                              !< If not provided, compute automatically
                num_min_kin=0                                                 !< Initialize kinetic mineral counter to zero
                num_sp=this%chem_syst%speciation_alg%num_aq_prim_species      !< Start offset after primary aqueous species
                
                !> Variable activity kinetic minerals - placed first after primary species
                do i=1,this%num_minerals_kin_var_act                          !< Loop through variable activity kinetic minerals
                    this%ind_min_Sk(num_min_kin+i)=num_sp+i                   !< Assign sequential index in kinetic matrix
                end do                                                        !< End loop
                num_min_kin=num_min_kin+this%num_minerals_kin_var_act         !< Update counter to reflect processed minerals
                
                !> Update offset to skip surface complexes, gases, secondary species
                num_sp=num_sp+this%num_minerals_kin_var_act+this%chem_syst%cat_exch_zone%num_surf_compl+&  !< Add variable act mins and surface complexes
                    this%chem_syst%gas_phase%num_var_act_species+this%chem_syst%speciation_alg%num_sec_aq_species  !< Add gases and secondary aqueous species
                
                !> Constant activity kinetic minerals - placed after variable activity
                do i=1,this%num_minerals_kin_cst_act                          !< Loop through constant activity kinetic minerals
                    this%ind_min_Sk(num_min_kin+i)=num_sp+this%ind_min_chem_syst(i)  !< Use index from chemical system with offset
                end do                                                        !< End loop
            end if                                                            !< End conditional block
        end subroutine
        
        !> \brief Set number of equilibrium minerals with variable activity
        !> \details
        !>   Sets count of equilibrium minerals whose activity varies with composition.
        !>   Variable activity equilibrium minerals contribute to mass action expressions.
        !> \param[inout] this Mineral zone object
        !> \param[in] num_mins_eq_var_act Count of equilibrium variable activity minerals
        subroutine set_num_mins_eq_var_act_min_zone(this,num_mins_eq_var_act)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: num_mins_eq_var_act                !< Equilibrium variable activity count
            this%num_minerals_eq_var_act=num_mins_eq_var_act                 !< Assign equilibrium variable activity count to member variable
        end subroutine
        
        !> \brief Set zone identifier
        !> \details
        !>   Sets unique identifier for this mineral zone.
        !>   Used to track zones during simulation.
        !> \param[inout] this Mineral zone object
        !> \param[in] id Zone identifier
        subroutine set_id(this,id)
            implicit none                                                     !< Disable implicit variable typing
            class(mineral_zone_c) :: this                                     !< Mineral zone object
            integer(kind=4), intent(in) :: id                                 !< Zone ID
            this%id=id                                                        !< Assign ID to member variable
        end subroutine
!> \brief End of mineral_zone_m module
end module