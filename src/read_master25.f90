!> @brief Reads 'master25_modif.dat' database
!> @details This subroutine parses the master25_modif.dat database file to populate
!> the chemical system object with species, minerals, gases, and surface complexes
!> @param[in,out] this Chemical system object to be populated
!> @param[in] path Directory path where the database file is located
!> @param[in] unit File unit number for reading the database
subroutine read_master25(this,path,unit)
    use chem_system_m, only: chem_system_c !< Import chemical system classes
    use array_ops_m, only: is_int_in_1D_array !< Import array operation utilities
    use aq_species_m, only: aq_species_c !< Import aqueous species class
    use gas_species_m, only: gas_species_c !< Import gas species class
    use mineral_m, only: mineral_c !< Import mineral class
    use solid_m, only: solid_species_c !< Import solid species class
    implicit none
    class(chem_system_c) :: this !< Chemical system object to be populated with database data
    character(len=*), intent(in) :: path !< Directory path containing the database file
    integer(kind=4), intent(in) :: unit !< File unit number for I/O operations
    
    !> @var wat_ind Index for water species
    !> @var i Generic loop counter
    !> @var j Generic loop counter for nested loops
    !> @var int_var Integer variable for I/O status checking
    !> @var num_reactants Number of reactants in current reaction
    !> @var ind_reacts Index for reactions array
    !> @var ind_mins Index for minerals array
    !> @var ind_gases Index for gases array
    !> @var ind_exch_cats Index for exchangeable cations array
    !> @var ind_cst_act_sp Index for constant activity species
    !> @var ind_var_act_sp Index for variable activity species
    !> @var num_temp_data Number of temperature data points in database
    !> @var ind_aq_compl Index for aqueous complexes
    !> @var num_pr_aq_sp Number of primary aqueous species (CHEPROO convention)
    !> @var num_aq_sp Total number of aqueous species
    !> @var num_aq_compl Number of aqueous complexes
    !> @var num_cst_act_sp Number of constant activity species
    !> @var num_var_act_sp Number of variable activity species
    !> @var num_eq_reacts Number of equilibrium reactions
    !> @var num_reacts Total number of reactions
    !> @var num_mins Number of minerals
    !> @var num_gases Number of gases
    !> @var num_exch_cats Number of exchangeable cations
    !> @var counter_aq_Compl Counter for aqueous complexes during second iteration
    !> @var counter_mins Counter for minerals during second iteration
    !> @var counter_surf_compl Counter for surface complexes during second iteration
    !> @var counter_gases Counter for gases during second iteration
    !> @var aq_sp_ind Index of aqueous species in species array
    !> @var min_ind Index of mineral in minerals array
    !> @var gas_ind Index of gas in gases array
    !> @var surf_compl_ind Index of surface complex in surface complexes array
    integer(kind=4) :: wat_ind,i,j,int_var,num_reactants,ind_reacts,ind_mins,ind_gases,ind_exch_cats,ind_cst_act_sp,&
        ind_var_act_sp,num_temp_data,ind_aq_compl,num_pr_aq_sp,num_aq_sp,num_aq_compl,num_cst_act_sp,num_var_act_sp,num_eq_reacts,&
        num_reacts,num_mins,num_gases,num_exch_cats,counter_aq_Compl,counter_mins,counter_surf_compl,counter_gases,aq_sp_ind,min_ind,&
        gas_ind,surf_compl_ind
    !> @var num_cst_act_gases Number of constant activity gases
    !> @var num_cst_act_mins Number of constant activity minerals
    !> @var num_var_act_gases Number of variable activity gases
    !> @var num_var_act_mins Number of variable activity minerals
    integer(kind=4) :: num_cst_act_gases,num_cst_act_mins,num_var_act_gases,num_var_act_mins
    !> @var react_indices Array of reaction indices for mapping database entries to reactions
    !> @var indices_exch_cats Array of indices for exchangeable cations
    !> @var indices_aq_sp Array of indices for aqueous species
    integer(kind=4), allocatable :: react_indices(:),indices_exch_cats(:),indices_aq_sp(:)
    !> @var sp_indices_chem_syst Temporary array to hold species indices in the chemical system for a reaction
    integer(kind=4), allocatable :: sp_indices_chem_syst(:)
    !> @var sp_ind_chem_syst Index of a species in the chemical system
    integer(kind=4) :: sp_ind_chem_syst
    !> @var aq_sp_flag Logical flag indicating if current species is an aqueous species
    !> @var min_flag Logical flag indicating if current species is a mineral
    !> @var gas_flag Logical flag indicating if current species is a gas
    !> @var surf_compl_flag Logical flag indicating if current species is a surface complex
    logical :: aq_sp_flag,min_flag,gas_flag,surf_compl_flag
    !> @var sp_flag Logical flag indicating if current species is found in the chemical system
    logical :: sp_flag
    !> @var log_K Logarithm (base 10) of equilibrium constant
    !> @var mol_weight Molecular weight in g/mol
    !> @var diff_vol Diffusion volume parameter
    !> @var mol_vol Molar volume in cm³/mol
    !> @var valence Valence (charge) of species
    real(kind=8) :: log_K,mol_weight,diff_vol,mol_vol,valence
    !> @var stoich_coeffs Array of stoichiometric coefficients for reaction
    real(kind=8), allocatable :: stoich_coeffs(:)
    !> @var str Temporary string variable for reading file
    !> @var name Name of current species being read
    !> @var filename Complete path and filename of database file
    character(len=256) :: str,name,filename
    !> @var species_names Array of species names participating in reaction
    character(len=256), allocatable :: species_names(:)
        
    !> @var aq_species Temporary aqueous species object for reading data
    type(aq_species_c) :: aq_species
    !> @var mineral Temporary mineral object for reading data
    type(mineral_c) :: mineral
    !> @var gas Temporary gas object for reading data
    type(gas_species_c) :: gas
    !> @var surf_compl Temporary surface complex object for reading data
    type(solid_species_c) :: surf_compl
    
    !> Construct full database filename by appending filename to path
    filename=trim(path)//'master25_modif.dat'

    !> Open database file for reading in old file mode (file must exist)
    open(unit,file=filename,status='old',action='read')
    
    !> Allocate arrays for tracking reaction indices and species indices during parsing
    allocate(react_indices(this%num_reacts),indices_aq_sp(this%aq_phase%num_species),indices_exch_cats(this%cat_exch_zone%num_exch_cats))
    
!> @section first_iteration First iteration of database
!> @details First pass through database to count species and allocate memory
    num_aq_sp=0 !< Initialize counter for total number of aqueous species
    num_pr_aq_sp=0 !< Initialize counter for number of primary aqueous species (CHEPROO)
    ind_aq_compl=0 !< Initialize index for aqueous complexes
    num_aq_compl=0 !< Initialize counter for number of aqueous complexes
    num_cst_act_sp=0 !< Initialize counter for number of constant activity species
    num_var_act_sp=0 !< Initialize counter for number of variable activity species
    num_cst_act_gases=0 !< Initialize counter for number of constant activity gases
    num_var_act_gases=0 !< Initialize counter for number of variable activity gases
    num_cst_act_mins=0 !< Initialize counter for number of constant activity minerals
    num_var_act_mins=0 !< Initialize counter for number of variable activity minerals
    num_eq_reacts=this%num_redox_eq_reacts !< Initialize counter for equilibrium reactions with redox reactions count
    num_reacts=this%num_redox_eq_reacts !< Initialize counter for total reactions with redox reactions count
    ind_mins=0 !< Initialize index for minerals
    num_mins=0 !< Initialize counter for number of minerals
    ind_gases=0 !< Initialize index for gases
    num_gases=0 !< Initialize counter for number of gases
    ind_exch_cats=0 !< Initialize index for exchangeable cations
    num_exch_cats=0 !< Initialize counter for number of exchangeable cations
    wat_ind=0 !< Initialize water index
    !> Read temperature line from database (string label and number of temperature data points)
    read(unit,*) str, num_temp_data
    
    !> @subsection primary_aq_species Primary aqueous species (CHEPROO)
    !> @details Loop through primary aqueous species in database
    do
        !> Read species data: name, ion size parameter, valence, molecular weight
        read(unit,*,iostat=int_var) aq_species%name, aq_species%params_act_coeff%ion_size_param, valence, mol_weight
        if (aq_species%name=='null') then !> Check if end of section marker is reached
            exit !> Exit loop when 'null' is encountered
        else
            !> Check if this species exists in the aqueous phase
            call this%aq_phase%is_species_in_aq_phase(aq_species,aq_sp_flag,aq_sp_ind)
            if (aq_sp_flag.eqv..true.) then !> If species is found in aqueous phase
                !> Set valence by converting real to integer
                call aq_species%set_valence(int(valence))
                !> Set molecular weight, converting from g/mol to kg/mol
                call aq_species%set_molecular_weight(mol_weight*1d-3)
                if (aq_species%name=='h2o(p)') then !> Check if species is water (primary)
                    !> Set constant activity flag to true for water
                    call aq_species%set_cst_act_flag(.true.)
                    !> Increment counter for constant activity species
                    num_cst_act_sp=num_cst_act_sp+1
                else !> For all other primary species
                    !> Set constant activity flag to false
                    call aq_species%set_cst_act_flag(.false.)
                    !> Legacy code for handling water index offset (commented out)
                    !if (aq_sp_ind>this%aq_phase%ind_wat .and. this%aq_phase%ind_wat>0) then !> chapuza
                    !    !aq_sp_ind=aq_sp_ind-1
                    !    call this%var_act_species(aq_sp_ind-1)%copy_species(aq_species)
                    !else
                    !    call this%var_act_species(aq_sp_ind)%copy_species(aq_species)
                    !end if
                    !> Increment counter for variable activity species
                    num_var_act_sp=num_var_act_sp+1
                end if
                !> Assign species to main species array
                call this%species(aq_sp_ind)%copy_species(aq_species)
                !> Store species index in tracking array
                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
                !> Assign species to aqueous phase species array
                call this%aq_phase%aq_species(aq_sp_ind)%copy_species(aq_species)
                !> Increment total aqueous species counter
                num_aq_sp=num_aq_sp+1
                !> Increment primary aqueous species counter
                num_pr_aq_sp=num_pr_aq_sp+1
            else !> If species not found in aqueous phase
                continue !> Skip to next iteration
            end if
        end if
    end do
    !> @subsection aqueous_complexes Aqueous complexes
    !> @details Loop through aqueous complexes to count and allocate reaction structures
    do
        !> Read complex name and number of reactants forming the complex
        read(unit,*,iostat=int_var) aq_species%name, num_reactants
        if (aq_species%name=='null') then !> Check if end of section marker is reached
            exit !> Exit loop when 'null' is encountered
        else
            !> Check if this complex exists in the aqueous phase
            call this%aq_phase%is_species_in_aq_phase(aq_species,aq_sp_flag,aq_sp_ind)
            if (aq_sp_flag.eqv..true.) then !> If complex is found in aqueous phase
                if (aq_species%name=='h2o') then !> Check if species is water
                    !> Set constant activity flag to true for water
                    call aq_species%set_cst_act_flag(.true.)
                    !> Increment constant activity species counter
                    num_cst_act_sp=num_cst_act_sp+1
                else !> For all other aqueous complexes
                    !> Set constant activity flag to false
                    call aq_species%set_cst_act_flag(.false.)
                    !> Increment variable activity species counter
                    num_var_act_sp=num_var_act_sp+1
                end if
            !> We assume all primary species are aqueous
                !> Allocate reaction with space for reactants plus product
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%allocate_reaction(num_reactants+1)
                !> Set reaction type to 1 (aqueous complexation)
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%set_react_type(1)
                !> Set reaction name to complex name
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%set_react_name(aq_species%name)
                !> Store reaction index for later data retrieval
                react_indices(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)=ind_aq_compl+1
                !> Store species index
                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
                !> Increment aqueous complexes counter
                num_aq_compl=num_aq_compl+1
                !> Increment total aqueous species counter
                num_aq_sp=num_aq_sp+1
            else !> If complex not found in aqueous phase
                continue !> Skip to next iteration
            end if
        end if
        !> Increment aqueous complex index
        ind_aq_compl=ind_aq_compl+1
    end do
    !> Update total equilibrium reactions count
    num_eq_reacts=num_eq_reacts+num_aq_compl
    !> Update total reactions count
    num_reacts=num_reacts+num_aq_compl
    !> @subsection minerals Minerals
    !> @details Loop through minerals to set properties and allocate reaction structures
    do
        !> Read mineral data: name, molar volume, and number of reactants
        read(unit,*,iostat=int_var) mineral%name, mol_vol, num_reactants
        if (mineral%name=='null') then !> Check if end of section marker is reached
            exit !> Exit loop when 'null' is encountered
        else
            !> Check if this mineral exists in the chemical system
            call this%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
            if (min_flag.eqv..true.) then !> If mineral is found in system
                !> Legacy code for assigning mineral species (commented out)
                !call mineral%mineral%copy_species(this%minerals(min_ind)%mineral)
                !call mineral%mineral%set_name(mineral%name)
                !> Set molar volume, converting from cm³/mol to m³/mol
                call this%minerals(min_ind)%mineral%set_mol_vol(mol_vol*1d-3)
                !> Legacy code for setting activity properties (commented out)
                !call mineral%mineral%set_cst_act_flag(.true.)
                !call mineral%mineral%set_valence(0)
                !this%minerals(min_ind)=mineral
                !> Assign mineral to main species array at appropriate index
                call this%species(this%aq_phase%num_species+min_ind)%copy_species(this%minerals(min_ind)%mineral)
                if (this%minerals(min_ind)%mineral%cst_act_flag.eqv..true.) then !> If mineral has constant activity
                    !> Increment constant activity minerals counter
                    num_cst_act_mins=num_cst_act_mins+1
                    !> Store species index in constant activity array
                    this%cst_act_sp_indices(num_cst_act_sp+num_cst_act_mins)=this%aq_phase%num_species+min_ind
                else !> If mineral has variable activity
                    !> Increment variable activity minerals counter
                    num_var_act_mins=num_var_act_mins+1
                    !> Store species index in variable activity array
                    this%var_act_sp_indices(num_var_act_sp+num_var_act_mins)=this%aq_phase%num_species+min_ind
                end if
                if (min_ind>this%num_minerals_kin) then !> If mineral undergoes equilibrium dissolution
                    !> Allocate equilibrium reaction with space for reactants plus product
                    call this%eq_reacts(num_eq_reacts+min_ind-this%num_minerals_kin)%allocate_reaction(num_reactants+1)
                    !> Set reaction type to 2 (mineral dissolution/precipitation)
                    call this%eq_reacts(num_eq_reacts+min_ind-this%num_minerals_kin)%set_react_type(2)
                    !> Set reaction name to mineral name
                    call this%eq_reacts(num_eq_reacts+min_ind-this%num_minerals_kin)%set_react_name(mineral%name)
                    !> Assign mineral species to last position in reaction species array
                    call this%eq_reacts(num_eq_reacts+min_ind-this%num_minerals_kin)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+min_ind)
                    !> Legacy code for setting stoichiometry (commented out)
                    !this%eq_reacts(num_eq_reacts+min_ind-this%num_minerals_kin)%stoichiometry(num_reactants+1)=-1d0
                else !> If mineral undergoes kinetic dissolution
                    !> Set kinetic reaction for this mineral
                    call this%kin_reacts(this%num_lin_kin_reacts+this%num_redox_kin_reacts+min_ind)%set_kin_reaction(&
                        this%min_kin_reacts(min_ind))
                    !> Allocate kinetic reaction with space for reactants plus product
                    call this%min_kin_reacts(min_ind)%allocate_reaction(num_reactants+1)
                    !> Set reaction type to 2 (mineral dissolution/precipitation)
                    call this%min_kin_reacts(min_ind)%set_react_type(2)
                    !> Set reaction name to mineral name
                    call this%min_kin_reacts(min_ind)%set_react_name(mineral%name)
                    !> Assign mineral species to last position in reaction species array
                    call this%min_kin_reacts(min_ind)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+min_ind)
                    !> Legacy code for setting stoichiometry (commented out)
                    !this%min_kin_reacts(min_ind)%stoichiometry(num_reactants+1)=-1d0
                end if
                !> Store reaction index for later data retrieval
                react_indices(num_reacts+min_ind)=ind_mins+1
                !> Increment minerals counter
                num_mins=num_mins+1
            else !> If mineral not found in system
                continue !> Skip to next iteration
            end if
        end if
        !> Increment mineral index
        ind_mins=ind_mins+1
    end do
    !> Update variable activity species count with variable activity minerals
    num_var_act_sp=num_var_act_sp+num_var_act_mins
    !> Update constant activity species count with constant activity minerals
    num_cst_act_sp=num_cst_act_sp+num_cst_act_mins
    !> Update total equilibrium reactions count with equilibrium minerals
    num_eq_reacts=num_eq_reacts+this%num_minerals_eq
    !> Update total reactions count with all minerals
    num_reacts=num_reacts+num_mins
    !> @subsection gases Gases
    !> @details Loop through gases to set properties and allocate reaction structures
    do
        !> Read gas data: name, molecular weight, diffusion volume, and number of reactants
        read(unit,*,iostat=int_var) gas%name, mol_weight, diff_vol, num_reactants
        if (gas%name=='null') then !> Check if end of section marker is reached
            exit !> Exit loop when 'null' is encountered
        else
            !> Check if this gas exists in the gas phase
            call this%gas_phase%is_gas_in_gas_phase(gas,gas_flag,gas_ind)
            if (gas_flag.eqv..true.) then !> If gas is found in gas phase
                !> Set molecular weight, converting from g/mol to kg/mol
                call this%gas_phase%gases(gas_ind)%set_molecular_weight(mol_weight*1d-3)
                !> TODO: Need to implement setter for diffusion volume
                !> falta set diff_vol
                !> Increment gases counter
                num_gases=num_gases+1
                !> Assign gas to main species array at appropriate index
                call this%species(this%aq_phase%num_species+this%num_minerals+gas_ind)%copy_species(gas)
                if (this%gas_phase%gases(gas_ind)%cst_act_flag.eqv..false.) then !> If gas has variable activity
                    !> Increment variable activity gases counter
                    num_var_act_gases=num_var_act_gases+1
                    !> Store species index in variable activity array
                    this%var_act_sp_indices(num_var_act_sp+num_var_act_gases)=this%aq_phase%num_species+this%num_minerals+gas_ind
                else !> If gas has constant activity
                    !> Increment constant activity gases counter
                    num_cst_act_gases=num_cst_act_gases+1
                    !> Store species index in constant activity array
                    this%cst_act_sp_indices(num_cst_act_sp+num_cst_act_gases)=this%aq_phase%num_species+this%num_minerals+gas_ind
                end if
                if (gas_ind<=this%gas_phase%num_gases_eq) then !> If gas undergoes equilibrium reaction
                    !> Allocate equilibrium reaction with space for reactants plus product
                    call this%eq_reacts(num_eq_reacts+gas_ind)%allocate_reaction(num_reactants+1)
                    !> Set reaction type to 6 (gas dissolution/exsolution)
                    call this%eq_reacts(num_eq_reacts+gas_ind)%set_react_type(6)
                    !> Set reaction name to gas name
                    call this%eq_reacts(num_eq_reacts+gas_ind)%set_react_name(gas%name)
                    !> Assign gas species to last position in reaction species array
                    call this%eq_reacts(num_eq_reacts+gas_ind)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+this%num_minerals+gas_ind)
                    !> Set stoichiometry to -1 for product
                    this%eq_reacts(num_eq_reacts+gas_ind)%stoichiometry(num_reactants+1)=-1d0
                    !> Store reaction index for later data retrieval
                    react_indices(num_reacts+gas_ind)=ind_gases+1
                end if
            else !> If gas not found in gas phase
                continue !> Skip to next iteration
            end if
        end if
        !> Increment gas index
        ind_gases=ind_gases+1
    end do
    !> Update variable activity species count with variable activity gases
    num_var_act_sp=num_var_act_sp+num_var_act_gases
    !> Update constant activity species count with constant activity gases
    num_cst_act_sp=num_cst_act_sp+num_cst_act_gases
    !> Update total equilibrium reactions count with gas reactions
    num_eq_reacts=num_eq_reacts+num_gases
    !> Update total reactions count with all gases
    num_reacts=num_reacts+num_gases
    !> @subsection surface_complexes Surface complexes
    !> @details Process surface complexes (exchangeable cations) if present in system
    if (this%cat_exch_zone%num_surf_compl>0) then !> Check if surface complexes exist
        !> Set name for exchange site species (X-)
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_name('x-')
        !> Set variable activity flag for exchange site
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_cst_act_flag(.false.)
        !> Set valence to -1 for exchange site
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_valence(-1)
        !> Store exchange site index in variable activity array
        this%var_act_sp_indices(num_var_act_sp+1)=this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1
        !> Legacy code for assigning exchange site to variable activity species (commented out)
        !call this%species(this%num_species-this%cat_exch%num_surf_compl+1)%copy_species(this%var_act_species(num_var_act_sp+1))
    end if
    !> Loop through surface complexes
    do
        !> Read surface complex name and number of reactants
        read(unit,*,iostat=int_var) surf_compl%name, num_reactants
        if (surf_compl%name=='null') then !> Check if end of section marker is reached
            exit !> Exit loop when 'null' is encountered
        else
            !> Check if this surface complex exists in cation exchange system
            call this%cat_exch_zone%is_surf_compl_in(surf_compl,surf_compl_flag,surf_compl_ind)
            if (surf_compl_flag.eqv..true.) then !> If surface complex is found
                !> Increment exchangeable cations counter
                num_exch_cats=num_exch_cats+1
                !> Legacy code for setting species properties (commented out)
                !call this%var_act_specis(num_var_act_sp+surf_compl_ind)%set_name(surf_compl%name)
                !call this%var_act_species(num_var_act_sp+surf_compl_ind)%set_cst_act_flag(.false.)
                !call this%var_act_species(num_var_act_sp+surf_compl_ind)%set_valence(0)
                !> Set name for surface complex in main species array
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_name(surf_compl%name)
                !> Set variable activity flag for surface complex
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_cst_act_flag(.false.)
                !> Set valence to 0 for surface complex
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_valence(0)
                !> Store species index in variable activity array
                this%var_act_sp_indices(num_var_act_sp+surf_compl_ind)=this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind
                !> Allocate equilibrium reaction with space for reactants plus product (offset by 1 for X-)
                call this%eq_reacts(num_eq_reacts+surf_compl_ind-1)%allocate_reaction(num_reactants+1)
                !> Set reaction type to 3 (surface complexation)
                call this%eq_reacts(num_eq_reacts+surf_compl_ind-1)%set_react_type(3)
                !> Store reaction index for later data retrieval (offset by 1 for X-)
                react_indices(num_reacts+surf_compl_ind-1)=ind_exch_cats+1
                !> Store surface complex index
                indices_exch_cats(num_exch_cats)=surf_compl_ind
            else !> If surface complex not found
                continue !> Skip to next iteration
            end if
        end if
        !> Increment exchangeable cations index
        ind_exch_cats=ind_exch_cats+1
    end do
    !> Update variable activity species count with surface complexes
    num_var_act_sp=num_var_act_sp+this%cat_exch_zone%num_surf_compl
    !> Update total equilibrium reactions count with surface complexation reactions
    num_eq_reacts=num_eq_reacts+num_exch_cats
    !> Update total reactions count with surface complexation reactions
    num_reacts=num_reacts+num_exch_cats
    
    !> Rewind file to beginning for second iteration
    rewind(unit)
    
!> @section second_iteration Second iteration of database
!> @details Second pass through database to read detailed reaction data (stoichiometry, equilibrium constants)
    counter_aq_compl=1 !< Initialize counter for aqueous complexes in second iteration
    ind_aq_compl=0 !< Reset index for aqueous complexes
    ind_cst_act_sp=0 !< Reset index for constant activity species
    ind_var_act_sp=num_pr_aq_sp-this%aq_phase%wat_flag !< Initialize index for variable activity species, accounting for water
    ind_mins=0 !< Reset index for minerals
    counter_mins=1 !< Initialize counter for minerals in second iteration
    ind_gases=0 !< Reset index for gases
    counter_gases=1 !< Initialize counter for gases in second iteration
    ind_exch_cats=0 !< Reset index for exchangeable cations
    counter_surf_compl=1 !< Initialize counter for surface complexes in second iteration
    !> Read and skip temperature line
    read(unit,*) str, num_temp_data
    !> @subsection primary_aq_species_skip Primary aqueous species
    !> @details Skip primary aqueous species section (already processed in first iteration)
    do
        !> Read species name to skip through section
        read(unit,*,iostat=int_var) name
        if (name=='null') exit !> Exit when end marker is reached
    end do
    !> @subsection aqueous_complexes_data Aqueous complexes
    !> @details Read detailed reaction data for aqueous complexes
    if (num_aq_compl>0) then !> Check if there are aqueous complexes to process
        do
            !> Check if current counter corresponds to a valid reaction index
            call is_int_in_1D_array(counter_aq_compl,react_indices(this%num_redox_eq_reacts+1:&
            this%num_redox_eq_reacts+num_aq_compl),aq_sp_flag,ind_reacts)
            if (aq_sp_flag.eqv..true.) then !> If this complex should be processed
                !> Increment aqueous complex index
                ind_aq_compl=ind_aq_compl+1
                !> Offset reaction index by number of redox equilibrium reactions
                ind_reacts=ind_reacts+this%num_redox_eq_reacts
                !> Allocate arrays for stoichiometry and species names
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read reaction data: complex name, number of reactants, reactant stoichiometry and names, log K, ion size, valence
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    num_reactants, (stoich_coeffs(j), species_names(j),j=1,num_reactants),&
                    log_K, aq_species%params_act_coeff%ion_size_param, valence!, mol_weight
                !> Set stoichiometric coefficients for the reaction
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                !> Set species names for the reaction
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                !> Set valence by converting real to integer
                call aq_species%set_valence(int(valence))
                !> Set molecular weight, converting from g/mol to kg/mol
                call aq_species%set_molecular_weight(mol_weight*1d-3)
                call aq_species%set_name(species_names(this%eq_reacts(ind_reacts)%num_species)) !> Set species name to reaction product name
                if (aq_species%name=='h2o') then !> Check if species is water
                    !> Increment constant activity species index
                    ind_cst_act_sp=ind_cst_act_sp+1
                    !> Set constant activity flag to true for water
                    call aq_species%set_cst_act_flag(.true.)
                    !> Legacy code for assigning to constant activity species array (commented out)
                    !call this%cst_act_species(1)%copy_species(aq_species)
                else !> For all other aqueous complexes
                    !> Increment variable activity species index
                    ind_var_act_sp=ind_var_act_sp+1
                    !> Set constant activity flag to false
                    call aq_species%set_cst_act_flag(.false.)
                    !> Legacy code for handling water index offset (commented out)
                    !if (indices_aq_sp(num_pr_aq_sp+ind_aq_compl)>this%aq_phase%ind_wat) then
                    !    call this%var_act_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl)-1)%copy_species(aq_species)
                    !else
                    !    call this%var_act_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%copy_species(aq_species)
                    !end if
                end if
                !> Assign species to main species array
                call this%species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%copy_species(aq_species)
                !> Assign species to aqueous phase species array
                call this%aq_phase%aq_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%copy_species(aq_species)
                !> Set equilibrium constant (converting from log K to K)
                call this%eq_reacts(ind_reacts)%set_eq_cst(10**(-log_K))
                !> Assign product species to last position in reaction species array
                call this%eq_reacts(ind_reacts)%set_single_species(&
                    this%eq_reacts(ind_reacts)%num_species, indices_aq_sp(num_pr_aq_sp+ind_aq_compl))
                !> Set stoichiometry to -1 for product
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                !> Change sign of all stoichiometry coefficients (convention)
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                !> Deallocate temporary arrays
                deallocate(stoich_coeffs,species_names)
            else !> If this entry should be skipped
                !> Read and discard line
                read(unit,*,iostat=int_var) name
                if (name=='null') exit !> Exit if end marker is reached
            end if
            !> Increment aqueous complex counter
            counter_aq_compl=counter_aq_compl+1
        end do
    else !> If no aqueous complexes exist
        !> Skip through section to reach end marker
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit !> Exit when end marker is reached
        end do
    end if
    !> @subsection minerals_data Minerals
    !> @details Read detailed reaction data for minerals
    if (num_mins>0) then !> Check if there are minerals to process
        do
            !> Check if current counter corresponds to a valid reaction index
            call is_int_in_1D_array(counter_mins,react_indices(this%num_redox_eq_reacts+num_aq_compl+1:this%num_redox_eq_reacts+ &
            num_aq_compl+num_mins),min_flag,ind_reacts)
            if ((min_flag .eqv. .true.) .AND. (ind_reacts>this%num_minerals_kin)) then !> If equilibrium mineral
                !> Offset reaction index to account for redox, complexes, and kinetic minerals
                ind_reacts=ind_reacts+this%num_redox_eq_reacts+num_aq_compl-this%num_minerals_kin
                !> Allocate arrays for stoichiometry and species names
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read reaction data: mineral name, molar volume, number of reactants, reactant stoichiometry and names, log K
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    mol_vol, num_reactants, (stoich_coeffs(j), species_names(j), j=1,num_reactants),&
                    log_K
                ! Legacy code for direct reading into reaction (commented out)
                ! read(unit,*,iostat=int_var) name, mol_vol, num_reactants, ((this%eq_reacts(ind_reacts)%stoichiometry(j), &

                ! this%eq_reacts(ind_reacts)%species(j)%name), j=1,num_reactants), log_K
                !> Set stoichiometric coefficients for the reaction
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                !> Set species names for the reaction
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(&
                    species_names,this%species)
                !> Set equilibrium constant (converting from log K to K)
                call this%eq_reacts(ind_reacts)%set_eq_cst(10**(-log_K))
                !> Set product mineral name on the species in the chemical system via species_ind
                sp_ind_chem_syst=this%eq_reacts(ind_reacts)%species_ind(this%eq_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Set stoichiometry to -1 for product
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                !> Change sign of all stoichiometry coefficients (convention)
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                !> Increment mineral index
                ind_mins=ind_mins+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) then !> If constant activity
                    !> Increment constant activity species index
                    ind_cst_act_sp=ind_cst_act_sp+1
                else !> If variable activity
                    !> Increment variable activity species index
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                !> Deallocate temporary arrays
                deallocate(stoich_coeffs,species_names)
            else if  (min_flag.eqv..true.) then !> If kinetic mineral
                !> Allocate arrays for stoichiometry and species names
                allocate(stoich_coeffs(this%min_kin_reacts(ind_reacts)%num_species),&
                    species_names(this%min_kin_reacts(ind_reacts)%num_species))
                !> Read reaction data: mineral name, molar volume, number of reactants, reactant stoichiometry and names, log K
                read(unit,*,iostat=int_var) species_names(this%min_kin_reacts(ind_reacts)%num_species), &
                    mol_vol, num_reactants, (stoich_coeffs(j), species_names(j), j=1,num_reactants),&
                    log_K
                ! Legacy code for direct reading into reaction (commented out)
                ! read(unit,*,iostat=int_var) name, mol_vol, num_reactants, ((this%min_kin_reacts(ind_reacts)%stoichiometry(j), &

                ! this%min_kin_reacts(ind_reacts)%species(j)%name), j=1,num_reactants), log_K
                !call this%kin_reacts(this%num_lin_kin_reacts+this%num_redox_kin_reacts+ind_reacts)%set_kin_reaction( &

                !    this%min_kin_reacts(ind_reacts))
                !> Set stoichiometric coefficients for the kinetic reaction
                call this%min_kin_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                !> Set species names for the kinetic reaction
                call this%min_kin_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                !> Set equilibrium constant (converting from log K to K)
                call this%min_kin_reacts(ind_reacts)%set_eq_cst(10**(-log_K))
                !> Set product mineral name on the species in the chemical system via species_ind
                sp_ind_chem_syst=this%min_kin_reacts(ind_reacts)%species_ind(this%min_kin_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%min_kin_reacts(ind_reacts)%num_species))
                !> Set stoichiometry to -1 for product
                this%min_kin_reacts(ind_reacts)%stoichiometry(this%min_kin_reacts(ind_reacts)%num_species)=-1d0
                !> Change sign of all stoichiometry coefficients (convention)
                call this%min_kin_reacts(ind_reacts)%change_sign_stoichiometry()
                !> Increment mineral index
                ind_mins=ind_mins+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) &
                then !> If constant activity
                    !> Increment constant activity species index
                    ind_cst_act_sp=ind_cst_act_sp+1
                else !> If variable activity
                    !> Increment variable activity species index
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                !> Deallocate temporary arrays
                deallocate(stoich_coeffs,species_names)
            else !> If this entry should be skipped
                !> Read and discard line
                read(unit,*,iostat=int_var) name
                if (name=='null') exit !> Exit if end marker is reached
            end if
            !> Increment minerals counter
            counter_mins=counter_mins+1
        end do
    else !> If no minerals exist
        !> Skip through section to reach end marker
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit !> Exit when end marker is reached
        end do
    end if
    !> @subsection gases_data Gases
    !> @details Read detailed reaction data for gases
    !> @note Gases should be ordered by activity for consistency with gas phase
    if (num_gases>0) then !> Check if there are gases to process
        do
            !> Check if current counter corresponds to a valid reaction index
            call is_int_in_1D_array(counter_gases,react_indices(this%num_redox_eq_reacts+num_aq_compl+num_mins+1:&
            this%num_redox_eq_reacts+num_aq_compl+num_mins+num_gases),gas_flag,ind_reacts)
            if (gas_flag.eqv..true.) then !> If this gas should be processed
                !> Offset reaction index to account for redox, complexes, and minerals
                ind_reacts=ind_reacts+this%num_redox_eq_reacts+num_aq_compl+this%num_minerals_eq
                !> Allocate arrays for stoichiometry and species names
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read reaction data: gas name, molecular weight, diffusion volume, number of reactants, reactant stoichiometry and names, log K
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    mol_weight, diff_vol, num_reactants, (stoich_coeffs(j), species_names(j),&
                    j=1,num_reactants), log_K
                ! Legacy code for direct reading into reaction (commented out)
                ! read(unit,*,iostat=int_var) name, mol_weight, diff_vol, num_reactants, &

                ! ((this%eq_reacts(ind_reacts)%stoichiometry(j), this%eq_reacts(ind_reacts)%species(j)%name), j=1,num_reactants), &

                ! log_K
                !> Set stoichiometric coefficients for the reaction
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                !> Set species names for the reaction
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                !> Set equilibrium constant (converting from log K to K)
                call this%eq_reacts(ind_reacts)%set_eq_cst(10**(-log_K))
                !> Set product gas properties on the species in the chemical system
                sp_ind_chem_syst=this%eq_reacts(ind_reacts)%species_ind(this%eq_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%eq_reacts(ind_reacts)%num_species))
                call this%species(sp_ind_chem_syst)%set_molecular_weight(mol_weight*1d-3)
                call this%species(sp_ind_chem_syst)%set_valence(0)
                !> Set stoichiometry to -1 for product
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                !> Change sign of all stoichiometry coefficients (convention)
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                !> Legacy code for incrementing variable activity index (commented out)
                !ind_var_act_sp=ind_var_act_sp+1
                !> Increment gas index
                ind_gases=ind_gases+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) then !> If constant activity
                    !> Increment constant activity species index
                    ind_cst_act_sp=ind_cst_act_sp+1
                else !> If variable activity
                    !> Increment variable activity species index
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                !> Deallocate temporary arrays
                deallocate(stoich_coeffs,species_names)
            else !> If this entry should be skipped
                !> Read and discard line
                read(unit,*,iostat=int_var) name
                if (name=='null') exit !> Exit if end marker is reached
            end if
            !> Increment gases counter
            counter_gases=counter_gases+1
        end do
    else !> If no gases exist
        !> Skip through section to reach end marker
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit !> Exit when end marker is reached
        end do
    end if
    !> @subsection surface_complexes_data Surface complexes
    !> @details Read detailed reaction data for surface complexes (exchangeable cations)
    if (num_exch_cats>0) then !> Check if there are surface complexes to process
        do
            !> Check if current counter corresponds to a valid reaction index
            call is_int_in_1D_array(counter_surf_compl,react_indices(num_aq_compl+num_mins+num_gases+1:num_aq_compl+num_mins+&
            num_gases+num_exch_cats),surf_compl_flag,ind_reacts)
            if (surf_compl_flag.eqv..true.) then !> If this surface complex should be processed
                !> Increment exchangeable cations index
                ind_exch_cats=ind_exch_cats+1
                !> Offset reaction index to account for redox, complexes, minerals, and gases
                ind_reacts=ind_reacts+this%num_redox_eq_reacts+num_aq_compl+this%num_minerals_eq+num_gases
                !> Allocate arrays for stoichiometry and species names
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read reaction data: surface complex name, number of reactants, reactant stoichiometry and names, log K, valence
                read(unit,*,iostat=int_var) surf_compl%name, num_reactants, (stoich_coeffs(j), species_names(j),j=1,num_reactants),&
                log_K, valence
                ! Legacy code for direct reading into reaction (commented out)
                ! read(unit,*,iostat=int_var) surf_compl%name, num_reactants, ((this%eq_reacts(ind_reacts)%stoichiometry(j),&

                !  this%eq_reacts(ind_reacts)%species(j)%name), j=1,num_reactants), log_K, surf_compl%valence
                !> Set stoichiometric coefficients for the reaction
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                !> Set species names for the reaction
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                !> Set valence by converting real to integer
                call surf_compl%set_valence(int(valence))
                !> Assign surface complex to cation exchange array
                this%cat_exch_zone%surf_compl(indices_exch_cats(ind_exch_cats))=surf_compl
                !> Assign surface complex to main species array at appropriate index
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+indices_exch_cats(ind_exch_cats))%copy_species(&
                    surf_compl)
                !> Legacy code for assigning to variable activity species (commented out)
                !call this%var_act_species(ind_var_act_sp+indices_exch_cats(ind_exch_cats))%copy_species(surf_compl)
                !> Set equilibrium constant (converting from log K to K)
                call this%eq_reacts(ind_reacts)%set_eq_cst(10**(-log_K))
                !> Set reaction name to surface complex name
                call this%eq_reacts(ind_reacts)%set_react_name(surf_compl%name)
                !> Assign product species to last position in reaction species array
                call this%eq_reacts(ind_reacts)%set_single_species(&
                    this%eq_reacts(ind_reacts)%num_species, &
                    this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+indices_exch_cats(ind_exch_cats))
                !> Set stoichiometry to -1 for product
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                !> Change sign of all stoichiometry coefficients (convention)
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                !> Deallocate temporary arrays
                deallocate(stoich_coeffs,species_names)
            else !> If this entry should be skipped
                !> Read and discard line
                read(unit,*,iostat=int_var) name
                if (name=='null') exit !> Exit if end marker is reached
            end if
            !> Increment surface complexes counter
            counter_surf_compl=counter_surf_compl+1
        end do
        !> Update variable activity species index with all surface complexes
        ind_var_act_sp=ind_var_act_sp+this%cat_exch_zone%num_surf_compl
    else !> If no surface complexes exist
        !> Skip through section to reach end marker
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit !> Exit when end marker is reached
        end do
    end if
    !> Close database file
    close(unit)
    !> Legacy code for rearranging species (commented out)
    !call this%rearrange_species()
    !> Legacy code for computing ionic strength squared (commented out)
    !call this%compute_z2() !> chapuza
end subroutine