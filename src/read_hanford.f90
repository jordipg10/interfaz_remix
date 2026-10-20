!> @file read_hanford.f90
!> @brief Reads the 'hanford.dat' multi-temperature thermodynamic database
!> @details This subroutine parses the hanford.dat database (EQ3/CHEPROO multi-temperature
!> format) to populate the chemical system object with species, aqueous complexes, gases,
!> minerals, and surface complexes. It is the multi-temperature analogue of read_master25.
!>
!> Key differences with respect to read_master25 (which reads master25_modif.dat):
!>   1. The header line lists a set of temperature points, e.g.
!>          'temperature points' 8 0. 25. 60. 100. 150. 200. 250. 300.
!>      and every reaction stores ONE log10(K) value PER temperature point (here 8).
!>      The reference-temperature value (25 C) is stored in each reaction's eq_cst_ref,
!>      and the FULL set of temperature-dependent constants is stored in the new
!>      reaction attribute eq_csts (see reaction_c). The reference index is located
!>      automatically as the temperature point closest to 25 C.
!>   2. The file section order is
!>          primary aqueous species -> aqueous complexes -> GASES -> MINERALS -> surface,
!>      i.e. gases come BEFORE minerals (the opposite of master25). The species/reaction
!>      ARRAY layout in the chemical system stays canonical (minerals then gases); this is
!>      achieved by computing the react_indices / eq_reacts base offsets from the structural
!>      counts (num_minerals, num_minerals_eq, gas_phase%num_species), so the processing order
!>      is irrelevant.
!>   3. The gas record uses a single leading field before the reactant count (like the
!>      minerals), and the molecular weight is the LAST field of the record (after the
!>      log10(K) block), for both gases and minerals.
!>
!> Equilibrium-constant sign convention follows read_master25: the database stores log10(K),
!> and REMIX stores K = 10**(-log10(K)) together with change_sign_stoichiometry(), so the
!> reaction is expressed with the product (the complex / mineral / gas / surface species) as
!> the species whose stoichiometric coefficient is -1.
!>
!> @param[in,out] this Chemical system object to be populated
!> @param[in]     path Directory path where the database file is located
!> @param[in]     unit File unit number for reading the database
!> @author Jordi Petchamé-Guerrero
subroutine read_hanford(this,path,unit)
    use chem_system_m, only: chem_system_c !< Import chemical system class
    use array_ops_m, only: is_int_in_1D_array !< Import array operation utilities
    use aq_species_m, only: aq_species_c !< Import aqueous species class
    use gas_species_m, only: gas_species_c !< Import gas species class
    use mineral_m, only: mineral_c !< Import mineral class
    use solid_m, only: solid_species_c !< Import solid species class
    implicit none
    class(chem_system_c) :: this !< Chemical system object to be populated with database data
    character(len=*), intent(in) :: path !< Directory path containing the database file
    integer(kind=4), intent(in) :: unit !< File unit number for I/O operations

    !> Loop counters, section indices and running counts (see read_master25 for details)
    integer(kind=4) :: wat_ind,i,j,int_var,num_reactants,ind_reacts,ind_mins,ind_gases,ind_exch_cats,ind_cst_act_sp,&
        ind_var_act_sp,num_temp_data,ind_aq_compl,num_pr_aq_sp,num_aq_sp,num_aq_compl,num_cst_act_sp,num_var_act_sp,num_eq_reacts,&
        num_reacts,num_mins,num_gases,num_exch_cats,counter_aq_Compl,counter_mins,counter_surf_compl,counter_gases,aq_sp_ind,min_ind,&
        gas_ind,surf_compl_ind
    integer(kind=4) :: num_cst_act_gases,num_cst_act_mins,num_var_act_gases,num_var_act_mins
    !> Multi-temperature bookkeeping
    integer(kind=4) :: ref_temp_ind !< Index (within the database temperature list) of the reference temperature (closest to 25 C)
    integer(kind=4) :: t !< Loop index over the database temperature points
    !> Canonical base offsets (computed from structural counts so the gas-before-mineral file
    !> order does not affect the array layout). "ri" = react_indices positions, "eq" = eq_reacts positions.
    integer(kind=4) :: base_ri_min,base_ri_gas,base_ri_surf !< react_indices bases for minerals / gases / surface
    integer(kind=4) :: base_eq_min,base_eq_gas,base_eq_surf !< eq_reacts bases for minerals / gases / surface
    integer(kind=4), allocatable :: react_indices(:),indices_exch_cats(:),indices_aq_sp(:)
    integer(kind=4), allocatable :: sp_indices_chem_syst(:)
    integer(kind=4) :: sp_ind_chem_syst
    logical :: aq_sp_flag,min_flag,gas_flag,surf_compl_flag
    logical :: sp_flag
    real(kind=8) :: log_K,mol_weight,diff_vol,mol_vol,valence
    real(kind=8), allocatable :: stoich_coeffs(:)
    real(kind=8), allocatable :: temp_points(:) !< Database temperature points [C]
    real(kind=8), allocatable :: logK_all(:) !< log10(K) at every temperature point for the current reaction
    real(kind=8), allocatable :: eq_csts_local(:) !< K = 10**(-log10(K)) at every temperature point for the current reaction
    character(len=256) :: str,name,filename
    character(len=256), allocatable :: species_names(:)

    type(aq_species_c) :: aq_species !< Temporary aqueous species object for reading data
    type(mineral_c) :: mineral !< Temporary mineral object for reading data
    type(gas_species_c) :: gas !< Temporary gas object for reading data
    type(solid_species_c) :: surf_compl !< Temporary surface complex object for reading data

    !> Construct full database filename by appending filename to path
    filename=trim(path)//'hanford.dat'

    !> Open database file for reading in old file mode (file must exist)
    open(unit,file=filename,status='old',action='read')

    !> Allocate arrays for tracking reaction indices and species indices during parsing
    allocate(react_indices(this%num_reacts),indices_aq_sp(this%aq_phase%num_species),indices_exch_cats(this%cat_exch_zone%num_exch_cats))

!> @section header Header: temperature points
!> @details Read the number of temperature points and the temperature list, then locate the
!> reference-temperature index (the point closest to 25 C). A backspace is used so the full
!> list can be read once the count is known.
    read(unit,*) str, num_temp_data
    if (num_temp_data<1) error stop "read_hanford: invalid number of temperature points"
    backspace(unit)
    allocate(temp_points(num_temp_data),logK_all(num_temp_data),eq_csts_local(num_temp_data))
    read(unit,*) str, num_temp_data, (temp_points(t), t=1,num_temp_data)
    ref_temp_ind=minloc(abs(temp_points-25d0),dim=1) !< reference temperature = closest point to 25 C

!> @section first_iteration First iteration of database
!> @details First pass to count species/reactions and allocate reaction structures.
    num_aq_sp=0
    num_pr_aq_sp=0
    ind_aq_compl=0
    num_aq_compl=0
    num_cst_act_sp=0
    num_var_act_sp=0
    num_cst_act_gases=0
    num_var_act_gases=0
    num_cst_act_mins=0
    num_var_act_mins=0
    num_eq_reacts=this%num_redox_eq_reacts
    num_reacts=this%num_redox_eq_reacts
    ind_mins=0
    num_mins=0
    ind_gases=0
    num_gases=0
    ind_exch_cats=0
    num_exch_cats=0
    wat_ind=0

    !> @subsection primary_aq_species Primary aqueous species
    do
        !> Read species data: name, ion size parameter, valence, molecular weight
        read(unit,*,iostat=int_var) aq_species%name, aq_species%params_act_coeff%ion_size_param, valence, mol_weight
        if (aq_species%name=='null') then !> End of section marker
            exit
        else
            call this%aq_phase%is_species_in_aq_phase(aq_species,aq_sp_flag,aq_sp_ind)
            if (aq_sp_flag.eqv..true.) then
                call aq_species%set_valence(int(valence))
                call aq_species%set_molecular_weight(mol_weight*1d-3) !> g/mol -> kg/mol
                if (aq_species%name=='h2o(p)' .or. aq_species%name=='H2O' .or. aq_species%name=='h2o') then !> water (constant activity)
                    call aq_species%set_cst_act_flag(.true.)
                    num_cst_act_sp=num_cst_act_sp+1
                else
                    call aq_species%set_cst_act_flag(.false.)
                    num_var_act_sp=num_var_act_sp+1
                end if
                call this%species(aq_sp_ind)%copy_species(aq_species)
                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
                call this%aq_phase%aq_species(aq_sp_ind)%copy_species(aq_species)
                num_aq_sp=num_aq_sp+1
                num_pr_aq_sp=num_pr_aq_sp+1
            else
                continue
            end if
        end if
    end do

    !> @subsection aqueous_complexes Aqueous complexes
    do
        !> Read complex name and number of reactants forming the complex
        read(unit,*,iostat=int_var) aq_species%name, num_reactants
        if (aq_species%name=='null') then
            exit
        else
            call this%aq_phase%is_species_in_aq_phase(aq_species,aq_sp_flag,aq_sp_ind)
            if (aq_sp_flag.eqv..true.) then
                if (aq_species%name=='h2o' .or. aq_species%name=='H2O') then !> water (constant activity)
                    call aq_species%set_cst_act_flag(.true.)
                    num_cst_act_sp=num_cst_act_sp+1
                else
                    call aq_species%set_cst_act_flag(.false.)
                    num_var_act_sp=num_var_act_sp+1
                end if
            !> We assume all primary species are aqueous
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%allocate_reaction(num_reactants+1)
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%set_react_type(1)
                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)%set_react_name(aq_species%name)
                react_indices(aq_sp_ind-num_pr_aq_sp+this%num_redox_eq_reacts)=ind_aq_compl+1
                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
                num_aq_compl=num_aq_compl+1
                num_aq_sp=num_aq_sp+1
            else
                continue
            end if
        end if
        ind_aq_compl=ind_aq_compl+1
    end do
    num_eq_reacts=num_eq_reacts+num_aq_compl
    num_reacts=num_reacts+num_aq_compl

    !> Now that num_aq_compl is known, compute the canonical base offsets from structural counts.
    !> These are order-independent, so gases (read before minerals) and minerals still land in
    !> their canonical positions in react_indices / eq_reacts.
    base_ri_min =this%num_redox_eq_reacts+num_aq_compl                 !< react_indices: minerals occupy [base_ri_min+1 .. +num_minerals]
    base_ri_gas =base_ri_min+this%num_minerals                        !< react_indices: gases occupy [base_ri_gas+1 .. +num_gases]
    base_ri_surf=base_ri_gas+this%gas_phase%num_species              !< react_indices: surface occupy [base_ri_surf+1 .. +num_exch_cats]
    base_eq_min =this%num_redox_eq_reacts+num_aq_compl                 !< eq_reacts: eq mineral idx = base_eq_min + (min_ind-num_minerals_kin)
    base_eq_gas =base_eq_min+this%num_minerals_eq                     !< eq_reacts: eq gas idx = base_eq_gas + gas_ind
    base_eq_surf=base_eq_gas+this%gas_phase%num_species              !< eq_reacts: surface idx = base_eq_surf + (surf_compl_ind-1)

    !> @subsection gases Gases (read BEFORE minerals in hanford.dat)
    do
        !> Read gas data: name, leading field (diffusion volume, unused), and number of reactants
        read(unit,*,iostat=int_var) gas%name, diff_vol, num_reactants
        if (gas%name=='null') then
            exit
        else
            call this%gas_phase%is_gas_in_gas_phase(gas,gas_flag,gas_ind)
            if (gas_flag.eqv..true.) then
                num_gases=num_gases+1
                call this%species(this%aq_phase%num_species+this%num_minerals+gas_ind)%copy_species(gas)
                if (this%gas_phase%gases(gas_ind)%cst_act_flag.eqv..false.) then
                    num_var_act_gases=num_var_act_gases+1
                    this%var_act_sp_indices(num_var_act_sp+num_var_act_gases)=this%aq_phase%num_species+this%num_minerals+gas_ind
                else
                    num_cst_act_gases=num_cst_act_gases+1
                    this%cst_act_sp_indices(num_cst_act_sp+num_cst_act_gases)=this%aq_phase%num_species+this%num_minerals+gas_ind
                end if
                if (gas_ind<=this%gas_phase%num_gases_eq) then !> equilibrium gas reaction
                    call this%eq_reacts(base_eq_gas+gas_ind)%allocate_reaction(num_reactants+1)
                    call this%eq_reacts(base_eq_gas+gas_ind)%set_react_type(6)
                    call this%eq_reacts(base_eq_gas+gas_ind)%set_react_name(gas%name)
                    call this%eq_reacts(base_eq_gas+gas_ind)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+this%num_minerals+gas_ind)
                    this%eq_reacts(base_eq_gas+gas_ind)%stoichiometry(num_reactants+1)=-1d0
                    react_indices(base_ri_gas+gas_ind)=ind_gases+1
                end if
            else
                continue
            end if
        end if
        ind_gases=ind_gases+1
    end do
    num_var_act_sp=num_var_act_sp+num_var_act_gases
    num_cst_act_sp=num_cst_act_sp+num_cst_act_gases
    num_eq_reacts=num_eq_reacts+num_gases
    num_reacts=num_reacts+num_gases

    !> @subsection minerals Minerals (read AFTER gases in hanford.dat)
    do
        !> Read mineral data: name, molar volume, and number of reactants
        read(unit,*,iostat=int_var) mineral%name, mol_vol, num_reactants
        if (mineral%name=='null') then
            exit
        else
            call this%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
            if (min_flag.eqv..true.) then
                call this%minerals(min_ind)%mineral%set_mol_vol(mol_vol*1d-3) !> cm3/mol (DB) -> dm3/mol
                call this%species(this%aq_phase%num_species+min_ind)%copy_species(this%minerals(min_ind)%mineral)
                if (this%minerals(min_ind)%mineral%cst_act_flag.eqv..true.) then
                    num_cst_act_mins=num_cst_act_mins+1
                    this%cst_act_sp_indices(num_cst_act_sp+num_cst_act_mins)=this%aq_phase%num_species+min_ind
                else
                    num_var_act_mins=num_var_act_mins+1
                    this%var_act_sp_indices(num_var_act_sp+num_var_act_mins)=this%aq_phase%num_species+min_ind
                end if
                if (min_ind>this%num_minerals_kin) then !> equilibrium mineral
                    call this%eq_reacts(base_eq_min+min_ind-this%num_minerals_kin)%allocate_reaction(num_reactants+1)
                    call this%eq_reacts(base_eq_min+min_ind-this%num_minerals_kin)%set_react_type(2)
                    call this%eq_reacts(base_eq_min+min_ind-this%num_minerals_kin)%set_react_name(mineral%name)
                    call this%eq_reacts(base_eq_min+min_ind-this%num_minerals_kin)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+min_ind)
                else !> kinetic mineral
                    call this%kin_reacts(this%num_lin_kin_reacts+this%num_redox_kin_reacts+min_ind)%set_kin_reaction(&
                        this%min_kin_reacts(min_ind))
                    call this%min_kin_reacts(min_ind)%allocate_reaction(num_reactants+1)
                    call this%min_kin_reacts(min_ind)%set_react_type(2)
                    call this%min_kin_reacts(min_ind)%set_react_name(mineral%name)
                    call this%min_kin_reacts(min_ind)%set_single_species(&
                        num_reactants+1, this%aq_phase%num_species+min_ind)
                end if
                react_indices(base_ri_min+min_ind)=ind_mins+1
                num_mins=num_mins+1
            else
                continue
            end if
        end if
        ind_mins=ind_mins+1
    end do
    num_var_act_sp=num_var_act_sp+num_var_act_mins
    num_cst_act_sp=num_cst_act_sp+num_cst_act_mins
    num_eq_reacts=num_eq_reacts+this%num_minerals_eq
    num_reacts=num_reacts+num_mins

    !> @subsection surface_complexes Surface complexes
    if (this%cat_exch_zone%num_surf_compl>0) then
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_name('x-')
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_cst_act_flag(.false.)
        call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1)%set_valence(-1)
        this%var_act_sp_indices(num_var_act_sp+1)=this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+1
    end if
    do
        read(unit,*,iostat=int_var) surf_compl%name, num_reactants
        if (surf_compl%name=='null') then
            exit
        else
            call this%cat_exch_zone%is_surf_compl_in(surf_compl,surf_compl_flag,surf_compl_ind)
            if (surf_compl_flag.eqv..true.) then
                num_exch_cats=num_exch_cats+1
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_name(surf_compl%name)
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_cst_act_flag(.false.)
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind)%set_valence(0)
                this%var_act_sp_indices(num_var_act_sp+surf_compl_ind)=this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+surf_compl_ind
                call this%eq_reacts(base_eq_surf+surf_compl_ind-1)%allocate_reaction(num_reactants+1)
                call this%eq_reacts(base_eq_surf+surf_compl_ind-1)%set_react_type(3)
                react_indices(base_ri_surf+surf_compl_ind-1)=ind_exch_cats+1
                indices_exch_cats(num_exch_cats)=surf_compl_ind
            else
                continue
            end if
        end if
        ind_exch_cats=ind_exch_cats+1
    end do
    num_var_act_sp=num_var_act_sp+this%cat_exch_zone%num_surf_compl
    num_eq_reacts=num_eq_reacts+num_exch_cats
    num_reacts=num_reacts+num_exch_cats

    !> Rewind file to beginning for second iteration
    rewind(unit)

!> @section second_iteration Second iteration of database
!> @details Second pass to read detailed reaction data (stoichiometry, equilibrium constants).
    counter_aq_compl=1
    ind_aq_compl=0
    ind_cst_act_sp=0
    ind_var_act_sp=num_pr_aq_sp-this%aq_phase%wat_flag
    ind_mins=0
    counter_mins=1
    ind_gases=0
    counter_gases=1
    ind_exch_cats=0
    counter_surf_compl=1
    !> Read and skip header temperature line (temperature list already parsed above)
    read(unit,*) str, num_temp_data

    !> @subsection primary_aq_species_skip Primary aqueous species (skip, already processed)
    do
        read(unit,*,iostat=int_var) name
        if (name=='null') exit
    end do

    !> @subsection aqueous_complexes_data Aqueous complexes
    if (num_aq_compl>0) then
        do
            call is_int_in_1D_array(counter_aq_compl,react_indices(this%num_redox_eq_reacts+1:&
            this%num_redox_eq_reacts+num_aq_compl),aq_sp_flag,ind_reacts)
            if (aq_sp_flag.eqv..true.) then
                ind_aq_compl=ind_aq_compl+1
                ind_reacts=ind_reacts+this%num_redox_eq_reacts
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read: complex name, num reactants, (stoich, species) pairs, num_temp_data logK values, ion size, valence, molecular weight
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    num_reactants, (stoich_coeffs(j), species_names(j),j=1,num_reactants),&
                    (logK_all(t), t=1,num_temp_data), aq_species%params_act_coeff%ion_size_param, valence, mol_weight
                log_K=logK_all(ref_temp_ind) !> reference-temperature log10(K)
                eq_csts_local=10d0**(-logK_all) !> K at every temperature point (sign convention as read_master25)
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                call aq_species%set_valence(int(valence))
                call aq_species%set_molecular_weight(mol_weight*1d-3)
                call aq_species%set_name(species_names(this%eq_reacts(ind_reacts)%num_species))
                if (aq_species%name=='h2o' .or. aq_species%name=='H2O') then
                    ind_cst_act_sp=ind_cst_act_sp+1
                    call aq_species%set_cst_act_flag(.true.)
                else
                    ind_var_act_sp=ind_var_act_sp+1
                    call aq_species%set_cst_act_flag(.false.)
                end if
                call this%species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%copy_species(aq_species)
                call this%aq_phase%aq_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%copy_species(aq_species)
                call this%eq_reacts(ind_reacts)%set_eq_cst_ref(10**(-log_K))
                call this%eq_reacts(ind_reacts)%set_eq_csts(eq_csts_local)
                call this%eq_reacts(ind_reacts)%set_single_species(&
                    this%eq_reacts(ind_reacts)%num_species, indices_aq_sp(num_pr_aq_sp+ind_aq_compl))
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                deallocate(stoich_coeffs,species_names)
            else
                read(unit,*,iostat=int_var) name
                if (name=='null') exit
            end if
            counter_aq_compl=counter_aq_compl+1
        end do
    else
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit
        end do
    end if

    !> @subsection gases_data Gases (read BEFORE minerals)
    if (num_gases>0) then
        do
            call is_int_in_1D_array(counter_gases,react_indices(base_ri_gas+1:base_ri_gas+this%gas_phase%num_species),gas_flag,ind_reacts)
            if (gas_flag.eqv..true.) then
                ind_reacts=ind_reacts+base_eq_gas
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read: gas name, leading field (diffusion volume, unused), num reactants, (stoich, species) pairs, num_temp_data logK values, molecular weight
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    diff_vol, num_reactants, (stoich_coeffs(j), species_names(j),j=1,num_reactants),&
                    (logK_all(t), t=1,num_temp_data), mol_weight
                log_K=logK_all(ref_temp_ind)
                eq_csts_local=10d0**(-logK_all)
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                call this%eq_reacts(ind_reacts)%set_eq_cst_ref(10**(-log_K))
                call this%eq_reacts(ind_reacts)%set_eq_csts(eq_csts_local)
                sp_ind_chem_syst=this%eq_reacts(ind_reacts)%species_ind(this%eq_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%eq_reacts(ind_reacts)%num_species))
                call this%species(sp_ind_chem_syst)%set_molecular_weight(mol_weight*1d-3)
                call this%species(sp_ind_chem_syst)%set_valence(0)
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                ind_gases=ind_gases+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) then
                    ind_cst_act_sp=ind_cst_act_sp+1
                else
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                deallocate(stoich_coeffs,species_names)
            else
                read(unit,*,iostat=int_var) name
                if (name=='null') exit
            end if
            counter_gases=counter_gases+1
        end do
    else
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit
        end do
    end if

    !> @subsection minerals_data Minerals (read AFTER gases)
    if (num_mins>0) then
        do
            call is_int_in_1D_array(counter_mins,react_indices(base_ri_min+1:base_ri_min+this%num_minerals),min_flag,ind_reacts)
            if ((min_flag .eqv. .true.) .AND. (ind_reacts>this%num_minerals_kin)) then !> equilibrium mineral
                min_ind=ind_reacts !> position of the mineral in this%minerals (kinetic minerals first)
                ind_reacts=ind_reacts+base_eq_min-this%num_minerals_kin
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read: mineral name, molar volume, num reactants, (stoich, species) pairs, num_temp_data logK values, molecular weight
                read(unit,*,iostat=int_var) species_names(this%eq_reacts(ind_reacts)%num_species), &
                    mol_vol, num_reactants, (stoich_coeffs(j), species_names(j), j=1,num_reactants),&
                    (logK_all(t), t=1,num_temp_data), mol_weight
                log_K=logK_all(ref_temp_ind)
                eq_csts_local=10d0**(-logK_all)
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                call this%eq_reacts(ind_reacts)%set_eq_cst_ref(10**(-log_K))
                call this%eq_reacts(ind_reacts)%set_eq_csts(eq_csts_local)
                sp_ind_chem_syst=this%eq_reacts(ind_reacts)%species_ind(this%eq_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%eq_reacts(ind_reacts)%num_species))
                !> molecular weight [kg/mol] on the mineral and on its chemical-system species (cumulative reaction amounts in kg)
                call this%minerals(min_ind)%mineral%set_molecular_weight(mol_weight*1d-3)
                call this%species(sp_ind_chem_syst)%set_molecular_weight(mol_weight*1d-3)
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                ind_mins=ind_mins+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) then
                    ind_cst_act_sp=ind_cst_act_sp+1
                else
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                deallocate(stoich_coeffs,species_names)
            else if (min_flag.eqv..true.) then !> kinetic mineral
                min_ind=ind_reacts !> kinetic minerals come first in this%minerals: the kinetic reaction index is the mineral's position
                allocate(stoich_coeffs(this%min_kin_reacts(ind_reacts)%num_species),&
                    species_names(this%min_kin_reacts(ind_reacts)%num_species))
                read(unit,*,iostat=int_var) species_names(this%min_kin_reacts(ind_reacts)%num_species), &
                    mol_vol, num_reactants, (stoich_coeffs(j), species_names(j), j=1,num_reactants),&
                    (logK_all(t), t=1,num_temp_data), mol_weight
                log_K=logK_all(ref_temp_ind)
                eq_csts_local=10d0**(-logK_all)
                call this%min_kin_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                call this%min_kin_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                call this%min_kin_reacts(ind_reacts)%set_eq_cst_ref(10**(-log_K))
                call this%min_kin_reacts(ind_reacts)%set_eq_csts(eq_csts_local)
                sp_ind_chem_syst=this%min_kin_reacts(ind_reacts)%species_ind(this%min_kin_reacts(ind_reacts)%num_species)
                call this%species(sp_ind_chem_syst)%set_name(species_names(this%min_kin_reacts(ind_reacts)%num_species))
                !> molecular weight [kg/mol] on the mineral and on its chemical-system species (cumulative reaction amounts in kg)
                call this%minerals(min_ind)%mineral%set_molecular_weight(mol_weight*1d-3)
                call this%species(sp_ind_chem_syst)%set_molecular_weight(mol_weight*1d-3)
                this%min_kin_reacts(ind_reacts)%stoichiometry(this%min_kin_reacts(ind_reacts)%num_species)=-1d0
                call this%min_kin_reacts(ind_reacts)%change_sign_stoichiometry()
                ind_mins=ind_mins+1
                if (this%species(sp_ind_chem_syst)%cst_act_flag.eqv..true.) then
                    ind_cst_act_sp=ind_cst_act_sp+1
                else
                    ind_var_act_sp=ind_var_act_sp+1
                end if
                deallocate(stoich_coeffs,species_names)
            else
                read(unit,*,iostat=int_var) name
                if (name=='null') exit
            end if
            counter_mins=counter_mins+1
        end do
    else
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit
        end do
    end if

    !> @subsection surface_complexes_data Surface complexes
    if (num_exch_cats>0) then
        do
            call is_int_in_1D_array(counter_surf_compl,react_indices(base_ri_surf+1:base_ri_surf+num_exch_cats),surf_compl_flag,ind_reacts)
            if (surf_compl_flag.eqv..true.) then
                ind_exch_cats=ind_exch_cats+1
                ind_reacts=ind_reacts+base_eq_surf
                allocate(stoich_coeffs(this%eq_reacts(ind_reacts)%num_species),&
                    species_names(this%eq_reacts(ind_reacts)%num_species))
                !> Read: surface complex name, num reactants, (stoich, species) pairs, num_temp_data logK values, valence
                read(unit,*,iostat=int_var) surf_compl%name, num_reactants, (stoich_coeffs(j), species_names(j),j=1,num_reactants),&
                (logK_all(t), t=1,num_temp_data), valence
                log_K=logK_all(ref_temp_ind)
                eq_csts_local=10d0**(-logK_all)
                species_names(this%eq_reacts(ind_reacts)%num_species)=surf_compl%name
                call this%eq_reacts(ind_reacts)%set_stoichiometry(stoich_coeffs)
                call this%eq_reacts(ind_reacts)%set_species_indices_from_names(species_names,this%species)
                call surf_compl%set_valence(int(valence))
                this%cat_exch_zone%surf_compl(indices_exch_cats(ind_exch_cats))=surf_compl
                call this%species(this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+indices_exch_cats(ind_exch_cats))%copy_species(&
                    surf_compl)
                call this%eq_reacts(ind_reacts)%set_eq_cst_ref(10**(-log_K))
                call this%eq_reacts(ind_reacts)%set_eq_csts(eq_csts_local)
                call this%eq_reacts(ind_reacts)%set_react_name(surf_compl%name)
                call this%eq_reacts(ind_reacts)%set_single_species(&
                    this%eq_reacts(ind_reacts)%num_species, &
                    this%speciation_alg%num_species-this%cat_exch_zone%num_surf_compl+indices_exch_cats(ind_exch_cats))
                this%eq_reacts(ind_reacts)%stoichiometry(this%eq_reacts(ind_reacts)%num_species)=-1d0
                call this%eq_reacts(ind_reacts)%change_sign_stoichiometry()
                deallocate(stoich_coeffs,species_names)
            else
                read(unit,*,iostat=int_var) name
                if (name=='null') exit
            end if
            counter_surf_compl=counter_surf_compl+1
        end do
        ind_var_act_sp=ind_var_act_sp+this%cat_exch_zone%num_surf_compl
    else
        do
            read(unit,*,iostat=int_var) name
            if (name=='null') exit
        end do
    end if

    !> Close database file
    close(unit)
end subroutine
