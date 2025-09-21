!> This class contains the chemical theory (like a "chemistry book"), that is, species, phases and reactions.
module chem_system_m
    use metodos_sist_lin_m
    use eq_reaction_m, only: eq_reaction_c, reaction_c
    use kin_reaction_m, only: kin_reaction_c, kin_reaction_poly_c, kin_params_c
    use lin_kin_reaction_m, only: lin_kin_reaction_c
    use redox_kin_reaction_m, only: redox_kin_c, Monod_params_c
    use kin_mineral_m, only: kin_mineral_c, kin_mineral_params_c, mineral_c
    use aq_phase_m, only: aq_phase_c, aq_species_c, phase_c
    use species_m, only: species_c, are_species_equal
    use gas_phase_m, only: gas_phase_c, gas_c
    use surf_compl_m, only: cat_exch_c, surface_c, solid_c
    use speciation_algebra_m, only: speciation_algebra_c
    use params_spec_vol_m, only: params_spec_vol_Redlich_c
    implicit none
    save
    type, public :: chem_system_c   !> chemical system class 
        integer(kind=4) :: num_species !> number of species
        integer(kind=4) :: num_reacts !> number of reactions
        type(species_c), allocatable :: species(:)  !> species
        real(kind=8), allocatable :: z2(:) !> squared charges of species
        type(aq_phase_c) :: aq_phase !> aqueous phase (we assume only 1)
        type(gas_phase_c) :: gas_phase !> gas phase (we assume only 1)
        integer(kind=4) :: num_var_act_species !> number of variable activity species
        integer(kind=4), allocatable :: var_act_sp_indices(:) !> variable activity species indices in "species" attribute
        !integer(kind=4) :: num_mineral_zones=0 !> number of mineral zones
        !type(mineral_zone_c), allocatable :: mineral_zones(:)
        integer(kind=4) :: num_minerals=0 !> number of minerals
        integer(kind=4) :: num_minerals_eq=0 !> number of minerals in equilibrium
        integer(kind=4) :: num_minerals_eq_cst_act=0 !> number of minerals in equilibrium with constant activity
        integer(kind=4) :: num_minerals_eq_var_act=0 !> number of minerals in equilibrium with variable activity
        integer(kind=4) :: num_minerals_kin=0 !> number of minerals NOT in equilibrium
        integer(kind=4) :: num_minerals_kin_cst_act=0 !> number of minerals NOT in equilibrium with constant activity
        integer(kind=4) :: num_minerals_kin_var_act=0 !> number of minerals NOT in equilibrium with variable activity
        type(mineral_c), allocatable :: minerals(:) !> minerals (first kinetic var act, then kin cst act, then eq var act, then equilibrium cst act)
        integer(kind=4) :: num_cst_act_species !> number of constant activity species
        integer(kind=4), allocatable :: cst_act_sp_indices(:) !> constant activity species indices in "species" attribute
        integer(kind=4) :: num_solids=0 !< number of solids
        type(cat_exch_c) :: cat_exch !> cation exchange object
        real(kind=8), allocatable :: stoich_mat(:,:) !> global stoichiometric matrix (S)
        real(kind=8), allocatable :: stoich_mat_sol(:,:) !> solid stoichiometric matrix (S_s)
        real(kind=8), allocatable :: stoich_mat_gas(:,:) !> gas stoichiometric matrix (S_g)
        real(kind=8), allocatable :: Se(:,:) !> equilibrium stoichiometric matrix
        real(kind=8), allocatable :: Sk(:,:) !> kinetic stoichiometric matrix
        integer(kind=4) :: num_eq_reacts=0 !> number of equilibrium reactions
        integer(kind=4) :: num_aq_eq_reacts=0 !> number of aqueous equilibrium reactions (we assume they are in equilibruim)
        integer(kind=4) :: num_aq_kin_reacts=0 !> number of aqueous kinetic reactions (we assume they are in equilibruim)
        integer(kind=4) :: num_redox_eq_reacts=0 !> number of redox equilibrium reactions
        type(eq_reaction_c), allocatable :: eq_reacts(:) !> equilibrium reactions
        type(speciation_algebra_c) :: speciation_alg !> speciation algebra object
        integer(kind=4) :: num_kin_reacts=0 !> number of kinetic reactions
        type(kin_reaction_poly_c), allocatable :: kin_reacts(:) !> kinetic reactions pointer array
        integer(kind=4) :: num_lin_kin_reacts=0 !> number of linear kinetic reactions
        type(lin_kin_reaction_c), allocatable :: lin_kin_reacts(:) !> linear kinetic reactions (we assume they are aqueous)
        !integer(kind=4) :: num_min_kin_reacts=0 !> number of mineral kinetic reactions
        type(kin_mineral_c), allocatable :: min_kin_reacts(:) !> mineral kinetic reactions
        integer(kind=4) :: num_redox_kin_reacts=0 !> number of Monod reactions
        type(redox_kin_c), allocatable :: redox_kin_reacts(:) !> redox kinetic reactions
        integer(kind=4) :: num_gas_kin_reacts=0 !> number of gas kinetic reactions
    contains
    !> Set
        procedure, public :: set_num_species
        procedure, public :: set_num_minerals
        procedure, public :: set_num_minerals_eq
        procedure, public :: set_num_minerals_eq_cst_act
        procedure, public :: set_num_minerals_eq_var_act
        procedure, public :: set_num_minerals_kin_cst_act
        procedure, public :: set_num_minerals_kin_var_act
        procedure, public :: set_num_eq_reacts
        procedure, public :: set_num_kin_reacts
        procedure, public :: set_num_gas_kin_reacts
        procedure, public :: set_num_lin_kin_reacts
        procedure, public :: set_num_minerals_kin
        procedure, public :: set_num_redox_kin_reacts
        procedure, public :: set_num_redox_eq_reacts
        procedure, public :: set_num_aq_eq_reacts
        procedure, public :: set_num_aq_kin_reacts
        procedure, public :: set_num_cst_act_species
        procedure, public :: set_num_var_act_species
        procedure, public :: set_species
        procedure, public :: set_cat_exch_obj
        procedure, public :: set_eq_reacts
        procedure, public :: set_kin_reacts
        procedure, public :: set_stoich_mat
        procedure, public :: set_stoich_mat_gas
        procedure, public :: set_stoich_mat_sol
    !> Allocate
        procedure, public :: allocate_species
        procedure, public :: allocate_cst_act_sp_indices
        procedure, public :: allocate_var_act_sp_indices
        procedure, public :: allocate_reacts
        procedure, public :: allocate_eq_reacts
        procedure, public :: allocate_kin_reacts
        procedure, public :: allocate_redox_kin_reacts
        procedure, public :: allocate_minerals
        procedure, public :: allocate_min_kin_reacts
        procedure, public :: allocate_lin_kin_reacts
    !> Compute
        procedure, public :: compute_num_kin_reacts
        procedure, public :: compute_num_reacts
        procedure, public :: compute_num_species
        procedure, public :: compute_z2
        procedure, public :: compute_num_solids_chem_syst
    !> Read
        !procedure, public :: read_chem_system
        procedure, public :: read_chem_system_CHEPROO
        procedure, public :: read_chem_system_PFLOTRAN
        procedure, public :: read_master25
        procedure, public :: read_kinetics_DB
        procedure, public :: read_Monod_reacts
        procedure, public :: read_PHREEQC_DB_opc1
        procedure, public :: read_PHREEQC_DB_opc2
    !> Is
        procedure, public :: is_mineral_in_chem_syst
        procedure, public :: is_eq_reaction_in_chem_syst
        procedure, public :: is_kin_reaction_in_chem_syst
        procedure, public :: is_species_in_chem_syst
        procedure, public :: are_minerals_equal
    !> Get
        procedure, public :: get_eq_csts
    !> Rearrange
        procedure, public :: rearrange_eq_reacts
        procedure, public :: rearrange_species
        !> Compute
        procedure, public :: compute_eq_csts_gases_cst_act
    end type
!****************************************************************************************************************************************************
    interface
        subroutine read_chem_system(this,path,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: filename
        end subroutine
        
        subroutine read_chem_system_CHEPROO(this,path_DB,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit
        end subroutine
        
        subroutine read_chem_system_PFLOTRAN(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
        end subroutine

        subroutine read_master25(this,path,unit)
            import chem_system_c
            import reaction_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
        end subroutine
        
        subroutine read_kinetics_DB(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
        end subroutine
        
        subroutine read_Monod_reacts(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: unit
        end subroutine
        
        subroutine read_PHREEQC_DB_opc1(this,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: filename
        end subroutine
        
        subroutine read_PHREEQC_DB_opc2(this,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: filename
        end subroutine
        
        subroutine read_PFLOTRAN_DB(this,unit,filename)
        import chem_system_c
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: filename
        end subroutine
        

        
      
       
     
       
    subroutine append_phase(this,phase)
        import chem_system_c
        import phase_c
        implicit none
        class(chem_system_c) :: this
        class(phase_c), intent(in) :: phase
    end subroutine
    
    subroutine append_reaction(this,reaction)
        import chem_system_c
        import reaction_c
        implicit none
        class(chem_system_c) :: this
        class(reaction_c), intent(in) :: reaction
    end subroutine
    
        subroutine set_stoich_mat(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
        end subroutine
        
        subroutine set_stoich_mat_gas(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
        end subroutine
        
        subroutine set_stoich_mat_sol(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
        end subroutine
      
    end interface
!***************************************************************************************************************************************************!
    contains
!*********************** SET ***********************************************************************************************************************!
        subroutine set_num_species(this,num_species)
        !< This subroutine sets the attribute "num_species"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                if (num_species<0) then
                    error stop "Number of species cannot be negative"
                else
                    this%num_species=num_species
                end if
            else
                this%num_species=size(this%species)
            end if
        end subroutine
        
        subroutine set_num_reacts(this,num_reacts)
        !< This subroutine sets the attribute "num_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_reacts
            if (present(num_reacts)) then
                if (num_reacts>this%num_species) then
                    error stop "Number of reactions cannot be greater than number of species"
                else if (num_reacts<0) then
                    error stop "Number of reactions cannot be negative"
                else
                    this%num_reacts=num_reacts
                end if
            else
                this%num_reacts=this%num_eq_reacts+this%num_kin_reacts
            end if
        end subroutine
        
        subroutine set_num_eq_reacts(this,num_eq_reacts)
        !< This subroutine sets the attribute "num_eq_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_eq_reacts
            this%num_eq_reacts=num_eq_reacts
        end subroutine

        subroutine set_num_kin_reacts(this,num_kin_reacts)
        !< This subroutine sets the attribute "num_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_kin_reacts
            if (num_kin_reacts>this%num_reacts) error stop "Number of kinetic reactions cannot be greater than number of reactions"
            this%num_kin_reacts=num_kin_reacts
        end subroutine

        subroutine set_num_gas_kin_reacts(this,num_gas_kin_reacts)
        !< This subroutine sets the attribute "num_gas_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_gas_kin_reacts
            if (num_gas_kin_reacts>this%num_kin_reacts) error stop "Number of gas kinetic reactions cannot be greater than number & 
                of kinetic reactions"
            this%num_gas_kin_reacts=num_gas_kin_reacts
        end subroutine
        
        subroutine set_species(this,species)
        !> This subroutine sets the "species" and "num_species" attributes 
            implicit none
            class(chem_system_c) :: this
            class(species_c), intent(in) :: species(:)
                        
            if (allocated(this%species) .and. size(species)/=this%num_species) then
                error stop "Wrong number of species"
            else
                this%species=species
                this%num_species=size(species)
            end if
        end subroutine
        
        subroutine set_num_cst_act_species(this,num_cst_act_species)
        !> This subroutine sets the "num_cst_act_species" attribute 
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_cst_act_species
            if (present(num_cst_act_species)) then
                this%num_cst_act_species=num_cst_act_species
            else
                this%num_cst_act_species=size(this%cst_act_sp_indices)
            end if
        end subroutine
        
        
        subroutine set_cat_exch_obj(this,cat_exch)
        !> This subroutine sets the "cat_exch" attribute 
            implicit none
            class(chem_system_c) :: this
            class(cat_exch_c), intent(in) :: cat_exch
            this%cat_exch=cat_exch
        end subroutine
                
        subroutine set_num_var_act_species(this,num_var_act_species)
        !> This subroutine sets the "num_var_act_species" attribute 
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_var_act_species
            if (present(num_var_act_species)) then
                this%num_var_act_species=num_var_act_species
            else if (allocated(this%var_act_sp_indices)) then
                this%num_var_act_species=size(this%var_act_sp_indices)
            else if (allocated(this%species) .and. allocated(this%cst_act_sp_indices)) then
                this%num_var_act_species=this%num_species-this%num_cst_act_species
            else
                error stop "Unable to compute the number of variable activity species"
            end if
        end subroutine
        
       
        
        subroutine set_num_minerals(this,num_minerals)
        !> This subroutine sets the "num_minerals"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals
            if (num_minerals<0) error stop "Number of minerals cannot be negative"
            this%num_minerals=num_minerals
        end subroutine
        
        subroutine set_num_minerals_eq(this,num_minerals_eq)
        !> This subroutine sets the "num_minerals_eq" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_eq
            if (num_minerals_eq>this%num_minerals .AND. allocated(this%minerals)) error stop "Number of minerals in equilibrium &
                cannot be greater than number of minerals"
            this%num_minerals_eq=num_minerals_eq
        end subroutine
        
        subroutine set_num_minerals_eq_cst_act(this,num_minerals_eq_cst_act)
        !> This subroutine sets the "num_minerals_eq_cst_act" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_eq_cst_act
            if (num_minerals_eq_cst_act>this%num_minerals_eq .AND. allocated(this%minerals)) then
                error stop "Number of minerals in equilibrium with constant activity cannot be greater than number of minerals in & 
                    equilibrium"
            end if
            this%num_minerals_eq_cst_act=num_minerals_eq_cst_act
        end subroutine
        
        subroutine set_num_minerals_eq_var_act(this,num_minerals_eq_var_act)
        !> This subroutine sets the "num_minerals_eq_var_act" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_eq_var_act
            if (num_minerals_eq_var_act>this%num_minerals_eq .AND. allocated(this%minerals)) then
                error stop "Number of minerals in equilibrium with constant activity cannot be greater than number of minerals in & 
                    equilibrium"
            end if
            this%num_minerals_eq_var_act=num_minerals_eq_var_act
        end subroutine
        
        subroutine set_num_minerals_kin_cst_act(this,num_minerals_kin_cst_act)
        !> This subroutine sets the "num_minerals_kin_cst_act" attribute
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_kin_cst_act
            if (num_minerals_kin_cst_act>this%num_minerals_kin .AND. allocated(this%minerals)) then
                error stop "Number of minerals NOT in equilibrium with constant activity cannot be greater than number of minerals &
                    NOT in equilibrium"
            end if
            this%num_minerals_kin_cst_act=num_minerals_kin_cst_act
        end subroutine
        
        subroutine set_num_minerals_kin_var_act(this,num_minerals_kin_var_act)
        !> This subroutine sets the "num_minerals_kin_var_act" attribute
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_kin_var_act
            if (num_minerals_kin_var_act>this%num_minerals_kin .AND. allocated(this%minerals)) then
                error stop "Number of minerals NOT in equilibrium with constant activity cannot be greater than number of minerals &
                    NOT in equilibrium"
            end if
            this%num_minerals_kin_var_act=num_minerals_kin_var_act
        end subroutine
        
        subroutine set_eq_reacts(this,eq_reacts)
        !> This subroutine sets the "eq_reacts" and "num_eq_reacts" attributes 
            implicit none
            class(chem_system_c) :: this
            class(eq_reaction_c), intent(in) :: eq_reacts(:)
            if (allocated(this%eq_reacts)) then
                deallocate(this%eq_reacts)
            end if
            this%eq_reacts=eq_reacts
            this%num_eq_reacts=size(this%eq_reacts)
        end subroutine
        
        subroutine set_kin_reacts(this,kin_reacts)
        !> This subroutine sets the "kin_reacts" and "num_kin_reacts" attributes 
            implicit none
            class(chem_system_c) :: this
            class(kin_reaction_poly_c), intent(in) :: kin_reacts(:)
            if (allocated(this%kin_reacts)) then
                deallocate(this%kin_reacts)
            end if
            this%kin_reacts=kin_reacts
            this%num_kin_reacts=size(this%kin_reacts)
        end subroutine
        
       subroutine set_num_lin_kin_reacts(this,num_lin_kin_reacts)
        !> This subroutine sets the "num_lin_kin_reacts" attribute 
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_lin_kin_reacts
            this%num_lin_kin_reacts=num_lin_kin_reacts
       end subroutine
       
       subroutine set_num_minerals_kin(this,num_minerals_kin)
        !> This subroutine sets the "num_lin_kin_reacts" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_minerals_kin
            this%num_minerals_kin=num_minerals_kin
       end subroutine
       
       subroutine set_num_redox_kin_reacts(this,num_redox_kin_reacts)
        !> This subroutine sets the "num_redox_kin_reacts" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_redox_kin_reacts
            this%num_redox_kin_reacts=num_redox_kin_reacts
       end subroutine
       
       subroutine set_num_redox_eq_reacts(this,num_redox_eq_reacts)
        !> This subroutine sets the "num_redox_eq_reacts" attribute
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_redox_eq_reacts
            this%num_redox_eq_reacts=num_redox_eq_reacts
       end subroutine
       
       subroutine set_num_aq_eq_reacts(this,num_aq_eq_reacts)
       !> This subroutine sets the "num_aq_eq_reacts" attribute
           class(chem_system_c) :: this
           integer(kind=4), intent(in), optional :: num_aq_eq_reacts
           if (present(num_aq_eq_reacts)) then
               this%num_aq_eq_reacts=num_aq_eq_reacts
           else
               this%num_aq_eq_reacts=this%num_redox_eq_reacts+this%aq_phase%num_aq_complexes
           end if
       end subroutine
       
       subroutine set_num_aq_kin_reacts(this,num_aq_kin_reacts)
       !> This subroutine sets the "num_aq_eq_reacts" attribute
           class(chem_system_c) :: this
           integer(kind=4), intent(in), optional :: num_aq_kin_reacts
           if (present(num_aq_kin_reacts)) then
               this%num_aq_kin_reacts=num_aq_kin_reacts
           else
               this%num_aq_kin_reacts=this%num_redox_kin_reacts+this%num_lin_kin_reacts
           end if
       end subroutine
       
!*********************** ALLOCATE ******************************************************************************************************************!
        subroutine allocate_species(this,num_species)
        !< This subroutine allocates the attribute "species"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                if (num_species<0) then
                    error stop "Number of species cannot be negative"
                else
                    this%num_species=num_species
                end if
            else
                call this%compute_num_species()
            end if
            if (allocated(this%species)) then
                deallocate(this%species)
            end if
            allocate(this%species(this%num_species))
        end subroutine
        
        subroutine allocate_cst_act_sp_indices(this,num_cst_act_species)
        !< This subroutine allocates the attribute "cst_act_sp_indices"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_cst_act_species
            if (present(num_cst_act_species)) then
                this%num_cst_act_species=num_cst_act_species                
            end if
            if (allocated(this%cst_act_sp_indices)) then
                deallocate(this%cst_act_sp_indices)
            end if
            allocate(this%cst_act_sp_indices(this%num_cst_act_species))
        end subroutine
        
        subroutine allocate_var_act_sp_indices(this,num_var_act_species)
        !< This subroutine allocates the attribute "var_act_sp_indices"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_var_act_species
            if (present(num_var_act_species)) then
                this%num_var_act_species=num_var_act_species                
            end if
            allocate(this%var_act_sp_indices(this%num_var_act_species))
        end subroutine
        
        subroutine allocate_reacts(this,num_eq_reacts,num_kin_reacts)
        !< This subroutine allocates the attributes "eq_reacts" & "kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in) :: num_eq_reacts,num_kin_reacts
            call this%allocate_eq_reacts(num_eq_reacts)
            call this%allocate_kin_reacts(num_kin_reacts)
            !call this%compute_num_reacts()
        end subroutine
        
        subroutine allocate_eq_reacts(this,num_eq_reacts)
        !< This subroutine allocates the attribute "eq_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_eq_reacts
            if (present(num_eq_reacts)) then
                this%num_eq_reacts=num_eq_reacts
            end if
            allocate(this%eq_reacts(this%num_eq_reacts))
        end subroutine
        
        subroutine allocate_kin_reacts(this,num_kin_reacts)
        !< This subroutine allocates the attribute "kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_kin_reacts
            if (present(num_kin_reacts)) then
                this%num_kin_reacts=num_kin_reacts
            end if
            allocate(this%kin_reacts(this%num_kin_reacts))
        end subroutine
        
        subroutine allocate_redox_kin_reacts(this,num_redox_kin_reacts)
        !< This subroutine allocates the attribute "redox_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_redox_kin_reacts
            if (present(num_redox_kin_reacts)) then
                this%num_redox_kin_reacts=num_redox_kin_reacts
            end if
            allocate(this%redox_kin_reacts(this%num_redox_kin_reacts))
        end subroutine
        
        subroutine allocate_lin_kin_reacts(this,num_lin_kin_reacts)
        !< This subroutine allocates the attribute "lin_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_lin_kin_reacts
            if (present(num_lin_kin_reacts)) then
                this%num_lin_kin_reacts=num_lin_kin_reacts
            end if
            allocate(this%lin_kin_reacts(this%num_lin_kin_reacts))
        end subroutine
        
        subroutine allocate_minerals(this,num_minerals)
        !< This subroutine allocates the attribute "minerals"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_minerals
            if (present(num_minerals)) then
                call this%set_num_minerals(num_minerals)
            end if
            allocate(this%minerals(this%num_minerals))
        end subroutine
        
        subroutine allocate_min_kin_reacts(this,num_min_kin)
        !< This subroutine allocates the attribute "min_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            integer(kind=4), intent(in), optional :: num_min_kin
            if (.not. allocated(this%minerals)) error stop
            if (present(num_min_kin)) then
                this%num_minerals_kin=num_min_kin
            end if
            allocate(this%min_kin_reacts(this%num_minerals_kin))
        end subroutine        

!*********************** COMPUTE ***********************************************************************************************************************!
        subroutine compute_num_species(this)
        !< This subroutine computes the attribute "num_species"
            implicit none
            class(chem_system_c) :: this
            this%num_species=this%num_var_act_species+this%num_cst_act_species
        end subroutine
        
        subroutine compute_num_kin_reacts(this)
        !< This subroutine computes the attribute "num_kin_reacts"
            implicit none
            class(chem_system_c) :: this
            this%num_kin_reacts=this%num_lin_kin_reacts+this%num_minerals_kin+this%num_redox_kin_reacts
        end subroutine
        
        
        
        
       
        
        
!********************************** IS *************************************************************************************************************!
        subroutine is_mineral_in_chem_syst(this,mineral,flag,mineral_ind)
        !> This subroutine checks if a mineral belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this                    !< chemical system
            class(mineral_c), intent(in) :: mineral                     !< mineral
            logical, intent(out) :: flag                                !< TRUE if mineral belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: mineral_ind       !> index of mineral in "minerals" attribute (if not belongs: 0)
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(mineral_ind)) then
                mineral_ind=0
            end if
            do i=1,this%num_minerals
                if (mineral%name==this%minerals(i)%name) then
                    flag=.true.
                    if (present(mineral_ind)) then
                        mineral_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
        subroutine is_species_in_chem_syst(this,species,flag,species_ind)
        !> This subroutine checks if a species belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            class(species_c), intent(in) :: species !< species
            logical, intent(out) :: flag !> TRUE if species belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: species_ind !> species index in "species" attribute (if not belongs: 0)
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(species_ind)) then
                species_ind=0
            end if
            do i=1,this%num_species
                if (species%name==this%species(i)%name) then
                    flag=.true.
                    if (present(species_ind)) then
                        species_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
        
        
        subroutine is_eq_reaction_in_chem_syst(this,react_name,flag,react_ind)
        !> This subroutine checks if an equilibrium reaction belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            character(len=*), intent(in) :: react_name !< reaction name
            logical, intent(out) :: flag !> TRUE if reaction belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: react_ind !< index in attribute "eq_reacts" (0 if not present)
            
            integer(kind=4) :: i,sp_ind
            integer(kind=4), allocatable :: sp_indices(:)
            logical :: sp_flag
            
            flag=.false.
            if (present(react_ind)) then
                react_ind=0
            end if
            do i=1,this%num_eq_reacts
                if (this%eq_reacts(i)%name==react_name) then
                    flag=.true.
                    if (present(react_ind)) then
                        react_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
        subroutine is_kin_reaction_in_chem_syst(this,react_name,flag,react_ind)
        !> This subroutine checks if a kinetic reaction belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            character(len=*), intent(in) :: react_name !< reaction name
            logical, intent(out) :: flag !> TRUE if reaction belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: react_ind !< index in attribute "kin_reacts" (0 if not present)
            
            integer(kind=4) :: i,sp_ind
            integer(kind=4), allocatable :: sp_indices(:)
            logical :: sp_flag
            
            flag=.false.
            if (present(react_ind)) then
                react_ind=0
            end if
            do i=1,this%num_kin_reacts
                if (this%kin_reacts(i)%kin_reaction%name==react_name) then
                    flag=.true.
                    if (present(react_ind)) then
                        react_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
!*********************** GET ***********************************************************************************************************************!
        function get_eq_csts(this) result(K)
        !> This function returns equilibrium constants of equilibrium reactions
            implicit none
            class(chem_system_c), intent(in) :: this
            real(kind=8), allocatable :: K(:)
            
            integer(kind=4) :: i
          
            allocate(K(this%num_eq_reacts))
            do i=1,this%num_eq_reacts
                K(i)=this%eq_reacts(i)%eq_cst
            end do
        end function
!*********************** REARRANGE *****************************************************************************************************************!
        subroutine rearrange_eq_reacts(this)
        !< This subroutine rearranges the "eq_reacts" attribute in the following order:
        !<      constant activity minerals in equilibrium
        !<      constant activity gases in equilibrium
        !<      redox equilibrium reactions
        !<      aqueous complexes
        !<      variable activity minerals in equilibrium
        !<      cation exchange
        !<      variable activity gases in equilibrium
            implicit none
            class(chem_system_c) :: this !> chemical system
            
            integer(kind=4) :: i,ind_min_cst_act,ind_aq,ind_gas_var_act,ind_surf,ind_gas_cst_act,ind_min_var_act,ind_redox
            type(eq_reaction_c), allocatable :: aux_eq_reacts(:)
            
            aux_eq_reacts=this%eq_reacts
            deallocate(this%eq_reacts)
            allocate(this%eq_reacts(this%num_eq_reacts))
        !> we initialise counters
            ind_min_cst_act=1                                                       !> constant activity minerals in equilibrium
            ind_gas_cst_act=ind_min_cst_act+this%num_minerals_eq_cst_act           !> constant activity gases in equilibrium
            ind_redox=ind_gas_cst_act+THIS%gas_phase%num_gases_eq_cst_act          !> redox equilibrium reactions
            ind_aq=ind_redox+this%num_redox_eq_reacts                                !> aqueous complexes
            ind_min_var_act=ind_aq+this%aq_phase%num_aq_complexes                            !> variable activity minerals
            ind_surf=ind_min_var_act+THIS%num_minerals_eq_var_act               !> cation exchange
            ind_gas_var_act=ind_surf+this%cat_exch%num_exch_cats                              !> variable activity gases
            do i=1,this%num_eq_reacts
                if (aux_eq_reacts(i)%react_type==2) then !> mineral dissolution/precipitation
                    if (aux_eq_reacts(i)%species(aux_eq_reacts(i)%num_species)%cst_act_flag.eqv..true.) then
                        this%eq_reacts(ind_min_cst_act)=aux_eq_reacts(i)
                        ind_min_cst_act=ind_min_cst_act+1
                    else
                        this%eq_reacts(ind_min_var_act)=aux_eq_reacts(i)
                        ind_min_var_act=ind_min_var_act+1
                    end if
                else if (aux_eq_reacts(i)%react_type==1) then !> aqueous complex
                    this%eq_reacts(ind_aq)=aux_eq_reacts(i)
                    ind_aq=ind_aq+1
                else if (aux_eq_reacts(i)%react_type==6) then !> gas
                    if (aux_eq_reacts(i)%species(aux_eq_reacts(i)%num_species)%cst_act_flag.eqv..true.) then
                        this%eq_reacts(ind_gas_cst_act)=aux_eq_reacts(i)
                        ind_gas_cst_act=ind_gas_cst_act+1
                    else
                        this%eq_reacts(ind_gas_var_act)=aux_eq_reacts(i)
                        ind_gas_var_act=ind_gas_var_act+1
                    end if
                else if (aux_eq_reacts(i)%react_type==3) then !> cation exchange
                    this%eq_reacts(ind_surf)=aux_eq_reacts(i)
                    ind_surf=ind_surf+1
                else if (aux_eq_reacts(i)%react_type==4) then !> redox
                    this%eq_reacts(ind_redox)=aux_eq_reacts(i)
                    ind_redox=ind_redox+1
                end if
            end do
        end subroutine
        
        subroutine rearrange_species(this)
        !< This subroutine rearranges the "species" attribute depending on the definition of the component matrix and the presence &
        !!    of surface complexes
            implicit none
            class(chem_system_c) :: this                                            !> chemical system
            
            integer(kind=4) :: i,num_sp,num_aq_sec,num_var_act_sp,num_cst_act_sp
            num_sp=0 !> counter number of species
            num_var_act_sp=0 !> counter number of variable activity species
            num_cst_act_sp=this%aq_phase%wat_flag !> counter number of constant activity species
            if ((this%speciation_alg%flag_comp .eqv. .false.) .and. (this%speciation_alg%flag_cat_exch .eqv. .true.)) then
            !<      primary aqueous species
            !<      free surface
            !<      var act minerals NOT in equilibrium
            !<      cst act minerals NOT in equilibrium
            !<      var act gases NOT in equilibrium
            !<      cst act gases NOT in equilibrium
            !<      aqueous complexes
            !<      var act minerals in equilibrium
            !<      cst act minerals in equilibrium
            !<      surface complexes
            !<      var act gases in equilibrium
            !<      cst act gases in equilibrium
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%assign_species(this%aq_phase%aq_species(i))
                end do
                call this%species(this%speciation_alg%num_aq_prim_species+1)%assign_species(this%cat_exch%surf_compl(1)) !> free surface
                num_sp=num_sp+this%speciation_alg%num_aq_prim_species+1
                do i=1,this%num_minerals_kin
                    call this%species(num_sp+i)%assign_species(this%minerals(i)%mineral)
                end do
                num_sp=num_sp+this%num_minerals_kin
                do i=1,this%gas_phase%num_gases_kin
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+i))
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin
                do i=1,this%aq_phase%num_aq_complexes
                    call this%species(num_sp+i)%assign_species(this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i))
                end do
                num_sp=num_sp+this%aq_phase%num_aq_complexes
                do i=1,this%num_minerals_eq
                    call this%species(num_sp+i)%assign_species(this%minerals(this%num_minerals_kin+i)%mineral)
                end do
                num_sp=num_sp+this%num_minerals_eq  
                do i=1,this%cat_exch%num_exch_cats
                    call this%species(num_sp+i)%assign_species(this%cat_exch%surf_compl(1+i))
                end do
                num_sp=num_sp+this%cat_exch%num_exch_cats
                do i=1,this%gas_phase%num_gases_eq
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(i))
                end do
            else if ((this%speciation_alg%flag_comp.eqv..true.) .and. (this%speciation_alg%flag_cat_exch.eqv..true.)) then
            !<      primary aqueous species
            !<      free surface
            !<      secondary variable activity aqueous species
            !<      kinetic variable activity minerals
            !<      equilibrium variable activity minerals
            !<      surface complexes
            !<      variable activity gases in equilibrium
            !<      ideal water
            !<      kinetic constant activity minerals
            !<      equilibrium constant activity minerals
            !<      constant activity gases in equilibrium
                call this%species(this%num_var_act_species+1)%assign_species(this%aq_phase%aq_species(this%aq_phase%ind_wat))
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%assign_species(this%aq_phase%aq_species(i))
                end do
                call this%species(this%speciation_alg%num_prim_species)%assign_species(this%cat_exch%surf_compl(1)) !> free surface
                do i=1,this%speciation_alg%num_aq_sec_var_act_species
                    call this%species(this%speciation_alg%num_prim_species+i)%assign_species(&
                    this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i))
                end do
                num_var_act_sp=1+this%speciation_alg%num_aq_var_act_species
                do i=1,this%num_minerals_kin_var_act
                    call this%species(num_var_act_sp+i)%assign_species(this%minerals(i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_kin_var_act
                do i=1,this%num_minerals_eq_var_act
                    call this%species(num_var_act_sp+i)%assign_species(this%minerals(this%num_minerals_kin+i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_eq_var_act
                ! do i=1,this%num_minerals
                !     if (this%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%assign_species(this%minerals(i)%mineral)
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%assign_species(this%minerals(i)%mineral)
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                do i=1,this%cat_exch%num_exch_cats
                    if (this%cat_exch%surf_compl(1+i)%cst_act_flag.eqv..false.) then
                        call this%species(num_var_act_sp+1)%assign_species(this%cat_exch%surf_compl(1+i))
                        num_var_act_sp=num_var_act_sp+1
                    else
                        call this%species(this%num_var_act_species+num_cst_act_sp+1)%assign_species(this%cat_exch%surf_compl(1+i))
                        num_cst_act_sp=num_cst_act_sp+1
                    end if
                end do
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%species(num_var_act_sp+i)%assign_species(this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                end do
                ! do i=1,this%gas_phase%num_species
                !     if (this%gas_phase%gases(i)%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%assign_species(this%gas_phase%gases(i))
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%assign_species(this%gas_phase%gases(i))
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                do i=1,this%num_minerals_kin_cst_act
                    call this%species(num_cst_act_sp+i)%assign_species(this%minerals(this%num_minerals_kin_var_act+i)%mineral)
                end do
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_kin_cst_act
                do i=1,this%num_minerals_eq_cst_act
                    call this%species(num_cst_act_sp+i)%assign_species(this%minerals(this%num_minerals_kin+&
                        this%num_minerals_eq_var_act+i)%mineral)
                end do
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_eq_cst_act
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    call this%species(num_cst_act_sp+i)%assign_species(this%gas_phase%gases(i))
                end do
            else if (this%speciation_alg%flag_comp .eqv. .true. .and. this%speciation_alg%flag_cat_exch .eqv. .false.) then
            !<      primary species
            !<      secondary variable activity aqueous species
            !<      kinetic variable activity minerals
            !<      equilibrium variable activity minerals
            !<      variable activity gases in equilibrium
            !<      ideal water
            !<      kinetic constant activity minerals
            !<      constant activity minerals in equilibrium
            !<      constant activity gases in equilibrium
                if (this%aq_phase%wat_flag==1) then
                    call this%species(this%num_var_act_species+1)%assign_species(this%aq_phase%aq_species(this%aq_phase%ind_wat))
                    this%cst_act_sp_indices(1)=this%num_var_act_species+1
                end if
                do i=1,this%speciation_alg%num_prim_species
                    call this%species(i)%assign_species(this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i)))
                    this%var_act_sp_indices(i)=i
                end do
                num_var_act_sp=num_var_act_sp+this%speciation_alg%num_prim_species
                do i=1,this%speciation_alg%num_aq_sec_var_act_species
                    call this%species(this%speciation_alg%num_prim_species+i)%assign_species(this%aq_phase%aq_species(&
                        this%aq_phase%ind_diss_solids(this%speciation_alg%num_aq_prim_species+i)))
                    this%var_act_sp_indices(num_var_act_sp+i)=this%speciation_alg%num_prim_species+i
                end do
                num_var_act_sp=num_var_act_sp+this%speciation_alg%num_aq_sec_var_act_species
                do i=1,this%num_minerals
                    if (this%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                        num_var_act_sp=num_var_act_sp+1
                        call this%species(num_var_act_sp)%assign_species(this%minerals(i)%mineral)
                        this%var_act_sp_indices(num_var_act_sp)=num_var_act_sp
                    else
                        num_cst_act_sp=num_cst_act_sp+1
                        call this%species(this%num_var_act_species+num_cst_act_sp)%assign_species(this%minerals(i)%mineral)
                        this%cst_act_sp_indices(num_cst_act_sp)=this%num_var_act_species+num_cst_act_sp
                    end if
                end do
                do i=1,this%gas_phase%num_species
                    if (this%gas_phase%gases(i)%cst_act_flag.eqv..false.) then
                        num_var_act_sp=num_var_act_sp+1
                        call this%species(num_var_act_sp)%assign_species(this%gas_phase%gases(i))
                        this%var_act_sp_indices(num_var_act_sp)=num_var_act_sp
                    else
                        num_cst_act_sp=num_cst_act_sp+1
                        call this%species(this%num_var_act_species+num_cst_act_sp)%assign_species(this%gas_phase%gases(i))
                        this%cst_act_sp_indices(num_cst_act_sp)=this%num_var_act_species+num_cst_act_sp
                    end if
                end do
            else
            !<      primary aqueous species
            !<      var act minerals NOT in equilibrium
            !<      cst act minerals NOT in equilibrium
            !<      gases var act NOT in equilibrium
            !<      gases cst act NOT in equilibrium
            !<      aqueous complexes
            !<      var act minerals in equilibrium
            !<      cst act minerals in equilibrium
            !<      var act gases in equilibrium
            !<      cst act gases in equilibrium
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%assign_species(this%aq_phase%aq_species(i))
                    if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .true.) then
                        this%cst_act_sp_indices(1)=i
                    else
                        num_var_act_sp=num_var_act_sp+1
                        this%var_act_sp_indices(num_var_act_sp)=i
                    end if
                end do
                num_sp=num_sp+this%speciation_alg%num_aq_prim_species
                do i=1,this%num_minerals_kin_var_act
                    call this%species(num_sp+i)%assign_species(this%minerals(i)%mineral)
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_kin_var_act
                num_sp=num_sp+this%num_minerals_kin_var_act
                do i=1,this%num_minerals_kin_cst_act
                    call this%species(num_sp+i)%assign_species(this%minerals(this%num_minerals_kin_var_act+i)%mineral)
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_kin_cst_act
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_kin_cst_act
                do i=1,this%gas_phase%num_gases_kin_var_act
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+&
                        this%gas_phase%num_gases_kin_cst_act+i))
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_kin_var_act
                num_sp=num_sp+this%gas_phase%num_gases_kin_var_act
                do i=1,this%gas_phase%num_gases_kin_cst_act
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+i))
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin_cst_act
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_kin_cst_act
                do i=1,this%speciation_alg%num_sec_aq_species
                    call this%species(num_sp+i)%assign_species(this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i))
                    if (this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i)%cst_act_flag .eqv. .true.) then
                        this%cst_act_sp_indices(1)=num_sp+i
                    else
                        num_var_act_sp=num_var_act_sp+1
                        this%var_act_sp_indices(num_var_act_sp)=num_sp+i
                    end if
                end do
                num_sp=num_sp+this%speciation_alg%num_sec_aq_species
                do i=1,this%num_minerals_eq_var_act
                    call this%species(num_sp+i)%assign_species(this%minerals(this%num_minerals_kin+i)%mineral)
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_eq_var_act
                num_var_act_sp=num_var_act_sp+this%num_minerals_eq_var_act
                do i=1,this%num_minerals_eq_cst_act
                    call this%species(num_sp+i)%assign_species(this%minerals(this%num_minerals_kin+&
                        this%num_minerals_eq_var_act+i)%mineral)
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_eq_cst_act
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_eq_cst_act
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_eq_var_act
                num_sp=num_sp+this%gas_phase%num_gases_eq_var_act
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    call this%species(num_sp+i)%assign_species(this%gas_phase%gases(i))
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%gas_phase%num_gases_eq_cst_act
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_eq_cst_act
            end if
        end subroutine
        
        subroutine compute_z2(this)
        !> This subroutine computes attribute "z2"
            implicit none
            class(chem_system_c) :: this !< chemical system
            
            integer(kind=4) :: i
            if (.not. allocated(this%z2)) then
                allocate(this%z2(this%num_species))
            end if
            do i=1,this%num_species
                this%z2(i)=this%species(i)%valence**2
            end do
        end subroutine
        
        subroutine compute_num_solids_chem_syst(this)
        !> This subroutine computes the "num_solids" attribute 
            implicit none
            class(chem_system_c) :: this
            this%num_solids=this%num_minerals+this%cat_exch%num_surf_compl
        end subroutine
        
        subroutine compute_num_reacts(this)
        !> This subroutine computes the "num_reacts" attribute 
            implicit none
            class(chem_system_c) :: this
            this%num_reacts=this%num_eq_reacts+this%num_kin_reacts
        end subroutine
        
        subroutine compute_eq_csts_gases_cst_act(this,part_pressures,ind_gases_eq_cst_act)
            class(chem_system_c) :: this
            real(kind=8), intent(in) :: part_pressures(:) !> partial pressures of gases in equilibrium with constant activity
            integer(kind=4), intent(in) :: ind_gases_eq_cst_act(:) !> indices of gases in equilibrium with constant activity in gas phase object
            
            integer(kind=4) :: i
            do i=1,size(ind_gases_eq_cst_act)
                this%eq_reacts(this%num_aq_eq_reacts+&
                    this%num_minerals_eq+ind_gases_eq_cst_act(i))%eq_cst = this%eq_reacts(&
                    this%num_aq_eq_reacts+this%num_minerals_eq + & 
                    ind_gases_eq_cst_act(i))%eq_cst/part_pressures(ind_gases_eq_cst_act(i))
            end do
        end subroutine
        
    subroutine are_minerals_equal(this,ind_mins_chem_syst_1,ind_mins_chem_syst_2,flag)
        !class(reactive_zone_c), intent(in) :: react_zone_1 !> first reactive zone
        !class(reactive_zone_c), intent(in) :: react_zone_2 !> second reactive zone
        class(chem_system_c), intent(in) :: this !> chemical system
        integer(kind=4), intent(in) :: ind_mins_chem_syst_1(:) !> indices of first set of minerals in chemical system
        integer(kind=4), intent(in) :: ind_mins_chem_syst_2(:) !> indices of second set of minerals in chemical system
        logical, intent(out) :: flag !> TRUE if minerals are equal, FALSE otherwise
        
        integer(kind=4) :: i
        logical :: min_flag !> mineral flag
 
        flag=.true. !> by default
        if (size(ind_mins_chem_syst_1)/=size(ind_mins_chem_syst_2)) then
            flag=.false.
        !else if (react_zone_1%num_cst_act_species/=react_zone_1%num_cst_act_species .or. react_zone_1%num_gases_eq_cst_act/=react_zone_1%num_gases_eq_cst_act) then
        !    flag=.false.
        else
            do i=1,size(ind_mins_chem_syst_1)
                call are_species_equal(this%minerals(ind_mins_chem_syst_1(i))%mineral,this%minerals(&
                    ind_mins_chem_syst_2(i))%mineral,min_flag)
                if (min_flag .eqv. .false.) then
                    flag=.false.
                    exit
                end if
            end do
        end if
    end subroutine

end module 