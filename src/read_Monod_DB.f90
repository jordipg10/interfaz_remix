!> @file read_Monod_DB.f90
!> @brief Reads Monod-type redox kinetic reactions from database file
!>
!> @details
!> This subroutine reads Monod reaction parameters from 'Monod_DB.dat' database file.
!> Monod reactions describe microbially-mediated redox reactions with Monod-type kinetics for
!> substrates and terminal electron acceptors (TEAs), including competitive inhibition.
!>
!> The Monod rate law implemented follows the form:
!> \f[
!> r = k \cdot \frac{[DOC]}{K_M^{DOC} + [DOC]} \cdot \frac{[TEA]}{K_M^{TEA} + [TEA]} \cdot \prod_{i=1}^{n_{inh}} \frac{K_i}{K_i + [I_i]}
!> \f]
!> where:
!> - \f$k\f$ = maximum rate constant
!> - \f$[DOC]\f$ = dissolved organic carbon concentration (electron donor)
!> - \f$K_M^{DOC}\f$ = half-saturation constant for DOC
!> - \f$[TEA]\f$ = terminal electron acceptor concentration
!> - \f$K_M^{TEA}\f$ = half-saturation constant for TEA
!> - \f$K_i\f$ = inhibition constant for inhibitor i
!> - \f$[I_i]\f$ = concentration of inhibitor i
!> - \f$n_{inh}\f$ = number of inhibitors
!>
!> File format:
!> ```
!> [First iteration - allocate arrays:]
!> reaction_label
!> n_inh
!> redox_couple_1_name  K_M_1
!> redox_couple_2_name  K_M_2
!> n_species
!> 
!> [Second iteration - read data:]
!> reaction_label
!> n_inh  inhibitor_1_name K_inh_1  inhibitor_2_name K_inh_2  ...
!> "stoichiometry"
!> n_species  species_1_name stoich_1  species_2_name stoich_2  ...  rate_cst
!> end
!> ```
!>
!> The subroutine performs two passes through the file:
!> 1. **First iteration:** Allocate arrays and read redox couple information
!> 2. **Second iteration:** Read full stoichiometry, inhibitors, and rate constants
!>
!> @param[in,out] this Chemical system object (modified with Monod reaction parameters)
!> @param[in] path Directory path where 'Monod_DB.dat' is located
!> @param[in] unit File unit number for reading
!>
!> @note DOC (dissolved organic carbon) species is assumed to be named 'ch2o(aq)'
!> @note Special handling for "aerobic degradation of DOC (review)" reaction with yield coefficients
!> @note The database file was created by jordi Petchamé-Guerrero
!> @note Reactions can be either kinetic reactions or equilibrium reactions
!>
!> @warning If inhibitor species not found in aqueous phase, program stops with error
!> @warning If redox couple species not found in aqueous phase, program stops with error
!> @warning Assumes inhibitors and redox couples are aqueous species (chapuza)
!>
!> @see redox_kin_c
!> @see Monod_params_c
!>
!> @author Jordi Petchamé-Guerrero
!> @date 2025
!>
!> Reads Monod reactions from database 'Monod_DB.dat' (created by me)
subroutine read_Monod_DB(this,path,unit)
    !> ================================================================
    !> Module imports - bring in required chemistry types
    !> ================================================================
    use chem_system_m, only: chem_system_c !> chemical system class
    use redox_kin_reaction_m, only: redox_kin_c
    use aq_species_m, only: aq_species_c
        !> redox_kin_c: redox kinetic reaction class
        !> Monod_params_c: Monod kinetic parameters class
        !> aq_species_c: aqueous species class
        !> kin_reaction_c: kinetic reaction base class
    implicit none !> enforce explicit variable declarations
    
    !> ================================================================
    !> Input/output parameters
    !> ================================================================
    class(chem_system_c) :: this !> chemical system object (modified - Monod reactions added)
    character(len=*), intent(in) :: path !> directory path to database file
    integer(kind=4), intent(in) :: unit !> file unit number for reading
            
    !> ================================================================
    !> Local variables - counters and indices
    !> ================================================================
    integer(kind=4) :: i,j,k,react_ind,n_sp,ind,num_inh,num_TEAs
        !> i: (not currently used) - general counter
        !> j: loop counter for species, inhibitors, redox couples
        !> k: (not currently used) - general counter
        !> react_ind: index of reaction in chemical system
        !> n_sp: number of species in reaction stoichiometry
        !> ind: index of species in aqueous phase
        !> num_inh: number of inhibitors (for equilibrium reactions)
        !> num_TEAs: number of terminal electron acceptors (for equilibrium reactions)
    
    !> Allocatable integer arrays
    integer(kind=4), allocatable :: n_inh(:),n_TEAs(:)
        !> n_inh: array storing number of inhibitors for each redox kinetic reaction
        !> n_TEAs: array storing number of terminal electron acceptors (not currently used)
    
    !> Real scalar variables
    real(kind=8) :: rate_cst,logK
        !> rate_cst: rate constant (for equilibrium reactions treated kinetically)
        !> logK: log10 of equilibrium constant
    
    !> Character variables for parsing
    character(len=256) :: label,filename,str
        !> label: reaction name/label read from file
        !> filename: full path to Monod reactions database file
        !> str: string for reading section headers (e.g., "stoichiometry")
    character(len=256), allocatable :: tmp_sp_names(:)
        !> tmp_sp_names: temporary array for reading species names from file
    
    !> Logical flag
    logical :: flag !> TRUE if reaction/species found in chemical system, FALSE otherwise
    
    !> ================================================================
    !> Temporary chemistry objects
    !> ================================================================
    type(aq_species_c) :: DOC !> dissolved organic carbon species (electron donor)
    type(aq_species_c), allocatable :: inhibitors(:) !> array of inhibitor species
    type(aq_species_c), allocatable :: redox_couple(:) !> array for redox couple (size 2: donor, acceptor)
    
!> ================================================================
!> SECTION 1: Pre-processing and initialization
!> ================================================================
!> Pre-process
    !> Construct full filename path
    filename=trim(path)//'\Monod_DB.dat' !> concatenate path + '\Monod_DB.dat'
    
    !> Allocate arrays
    allocate(n_inh(this%num_redox_kin_reacts)) !> allocate inhibitor count array (size = number of redox kinetic reactions)
    allocate(redox_couple(2)) !> allocate redox couple array (size 2: electron donor and acceptor)
    
    !> Set DOC species name
    call DOC%set_name('ch2o(aq)') !> DOC (dissolved organic carbon) is named 'ch2o(aq)'
    
!> ================================================================
!> SECTION 2: Main processing - two-pass file reading
!> ================================================================
!> Process
    !> Open Monod reactions database file for reading
    open(unit,file=filename,status='old',action='read') !> status='old' means file must exist
    
!> ================================================================
!> SECTION 3: First iteration - allocate arrays and read structure
!> ================================================================
!> First iteration: Count species, allocate arrays, read redox couples and half-saturation constants
    !react_ind=1 !> (COMMENTED) counter for redox kinetic reactions in chemical system
    do !> infinite loop - exit when 'end' encountered
        !> Read reaction label/name
        read(unit,*) label !> reaction name/identifier
        
        if (label=='end') then !> if end-of-data sentinel found
            rewind(unit) !> rewind file to beginning for second iteration
            exit !> exit first iteration loop
        else !> continue processing reaction
            !> ----------------------------------------------------------------
            !> Check if this is a redox kinetic reaction
            !> ----------------------------------------------------------------
            call this%is_kin_reaction_in_chem_syst(label,flag,react_ind)
                !> Search for reaction in kinetic reactions list
                !> Returns: flag (TRUE if found), react_ind (index in system)
            
            if (flag.eqv. .true.) then !> if reaction is a kinetic reaction
                react_ind=react_ind-this%num_lin_kin_reacts !> adjust index to redox kinetic reactions only
                !> ------------------------------------------------------------
                !> Read and setup kinetic reaction structure
                !> ------------------------------------------------------------
                read(unit,*) n_inh(react_ind) !> read number of inhibitors for this reaction
                
                !> Allocate inhibitor constants array
                call this%redox_kin_reacts(react_ind)%params%allocate_k_inh(n_inh(react_ind))
                    !> allocate k_inh array (size = n_inh)
                
                !> Allocate half-saturation constants array
                call this%redox_kin_reacts(react_ind)%params%allocate_k_M()
                    !> allocate k_M array (size 2: for DOC and TEA)
                
                !> Compute total number of terms (inhibitors + redox couples)
                call this%redox_kin_reacts(react_ind)%params%compute_num_terms()
                    !> num_terms = n_inh + 2 (two redox couples: DOC and TEA)
                
                !> Allocate indices array for aqueous phase species
                call this%redox_kin_reacts(react_ind)%allocate_indices_aq_phase_kin_react(&
                    this%redox_kin_reacts(react_ind)%params%num_terms)
                    !> allocate indices_aq_phase array (size = num_terms)
                !allocate(this%redox_kin_reacts(react_ind)%indices_aq_phase(this%redox_kin_reacts(react_ind)%params%num_terms)) !> (COMMENTED) alternative allocation
                
                !> Read redox couple names and half-saturation constants
                read(unit,*) (redox_couple(j)%name, this%redox_kin_reacts(react_ind)%params%k_M(j), j=1,2)
                    !> j=1: electron donor (DOC) and its K_M
                    !> j=2: terminal electron acceptor (TEA) and its K_M
                
                !> Find redox couple species indices in aqueous phase
                do j=1,2 !> loop over redox couple (donor and acceptor)
                    if (trim(redox_couple(j)%name)=='*') then !> placeholder: no redox couple
                        this%redox_kin_reacts(react_ind)%indices_aq_phase(this%redox_kin_reacts(react_ind)%params%n_inh+j)=0
                        cycle
                    end if
                    call this%aq_phase%is_species_in_aq_phase(redox_couple(j),flag,ind)
                        !> search for species in aqueous phase - chapuza (puede haber especies no acuosas)
                        !> Returns: flag (TRUE if found), ind (index in aqueous phase)
                    
                    if (flag.eqv..false.) then !> if redox couple species not found
                        error stop 'Redox couple not found in chemical system' !> stop with error
                    else !> species found
                        !> Store species index (after inhibitors in indices array)
                        this%redox_kin_reacts(react_ind)%indices_aq_phase(this%redox_kin_reacts(react_ind)%params%n_inh+j)=ind
                            !> indices layout: [inhibitors (1..n_inh) | DOC (n_inh+1) | TEA (n_inh+2)]
                    end if
                end do
                
                !> Read number of species in reaction stoichiometry
                read(unit,*) n_sp !> total number of species in reaction
                
                !> Allocate reaction arrays (species, stoichiometry)
                call this%redox_kin_reacts(react_ind)%allocate_reaction(n_sp)
                    !> allocate species array and stoichiometry array (size n_sp)
                
            else !> not a kinetic reaction - check if equilibrium reaction
                !> ----------------------------------------------------------------
                !> Check if this is an equilibrium reaction
                !> ----------------------------------------------------------------
                call this%is_eq_reaction_in_chem_syst(label,flag,react_ind)
                    !> Search for reaction in equilibrium reactions list
                
                if (flag .eqv. .true.) then !> if reaction is an equilibrium reaction
                    !> Read equilibrium reaction structure info
                    read(unit,*) num_inh !> read number of inhibitors (for equilibrium reaction)
                    read(unit,*) num_TEAs !> read number of terminal electron acceptors
                    read(unit,*) n_sp !> read number of species in stoichiometry
                    
                    !> Allocate equilibrium reaction arrays
                    call this%eq_reacts(react_ind)%allocate_reaction(n_sp)
                        !> allocate species and stoichiometry arrays
                else !> reaction not found in either kinetic or equilibrium reactions
                    continue !> skip this reaction, continue to next
                end if
            end if
        end if
    end do !> end first iteration loop
    
!> ================================================================
!> SECTION 4: Second iteration - read detailed reaction data
!> ================================================================
!> Second iteration: Read complete reaction data including inhibitors, species, stoichiometry, and rate constants
    !ind_redox_kin=1 !> (COMMENTED) counter for Monod reactions in chemical system
    do !> infinite loop - exit when 'end' encountered
        !> Read reaction label/name
        read(unit,*) label !> reaction identifier
        
        if (label=='end') then !> if end-of-data sentinel found
            exit !> exit second iteration loop - processing complete
        else !> continue processing reaction
            !> ----------------------------------------------------------------
            !> Check if this is a redox kinetic reaction
            !> ----------------------------------------------------------------
            call this%is_kin_reaction_in_chem_syst(label,flag,react_ind)
                !> Search for reaction in kinetic reactions list
                !> Returns: flag (TRUE if found), react_ind (index in chemical system)
            
            if (flag.eqv. .true.) then !> if reaction is a kinetic reaction
                react_ind=react_ind-this%num_lin_kin_reacts !> adjust index to redox kinetic reactions only
                !> ------------------------------------------------------------
                !> Read inhibitor data (conditionally based on inhibitor count)
                !> ------------------------------------------------------------
                if (n_inh(react_ind)>0) then !> if inhibitors present for this reaction
                    !> Allocate temporary inhibitor array
                    allocate(inhibitors(n_inh(react_ind))) !> allocate temporary storage for inhibitor species
                    
                    !> Read number of inhibitors, species names, and inhibition constants
                    read(unit,*) this%redox_kin_reacts(react_ind)%params%n_inh,&
                        (inhibitors(j)%name, this%redox_kin_reacts(react_ind)%params%k_inh(j), &
                        j=1,this%redox_kin_reacts(react_ind)%params%n_inh)
                        !> reads: n_inh, name_1, K_inh_1, name_2, K_inh_2, ..., name_n, K_inh_n
                    
                    !> Find inhibitor species indices in aqueous phase
                    do j=1,this%redox_kin_reacts(react_ind)%params%n_inh !> loop over all inhibitors
                        call this%aq_phase%is_species_in_aq_phase(inhibitors(j),flag,ind)
                            !> search for inhibitor species in aqueous phase - chapuza (puede haber especies no acuosas)
                            !> Returns: flag (TRUE if found), ind (index in aqueous phase)
                        
                        if (flag.eqv. .false.) then !> if inhibitor species not found in aqueous phase
                            error stop 'Inhibitor not found in chemical system' !> stop with error message
                        else !> inhibitor species found
                            !> Store inhibitor index in first positions of indices array
                            this%redox_kin_reacts(react_ind)%indices_aq_phase(j)=ind
                                !> indices layout: [inhibitors (1..n_inh) | DOC | TEA]
                        end if
                    end do
                    
                    !> Free temporary storage
                    deallocate(inhibitors) !> deallocate temporary inhibitor array
                else !> no inhibitors for this reaction
                    !> Read number of inhibitors (should be 0)
                    read(unit,*) this%redox_kin_reacts(react_ind)%params%n_inh
                        !> reads just n_inh value (expected to be 0)
                end if
                
                !> ------------------------------------------------------------
                !> Skip redox couple line (already read in first iteration)
                !> ------------------------------------------------------------
                read(unit,*) str !> read and discard redox couple line (already processed)
                
                !> ------------------------------------------------------------
                !> Read reaction stoichiometry, species, and rate constant
                !> ------------------------------------------------------------
                if (allocated(tmp_sp_names)) deallocate(tmp_sp_names)
                allocate(tmp_sp_names(this%redox_kin_reacts(react_ind)%num_species))
                read(unit,*) this%redox_kin_reacts(react_ind)%num_species, (tmp_sp_names(j), &
                this%redox_kin_reacts(react_ind)%stoichiometry(j), j=1,this%redox_kin_reacts(react_ind)%num_species),&
                this%redox_kin_reacts(react_ind)%params%rate_cst!, this%redox_kin_reacts(react_ind)%eq_cst
                    !> reads: n_sp, name_1, coeff_1, name_2, coeff_2, ..., name_n, coeff_n, k
                    !> coeff: negative=reactant, positive=product
                    !> k: kinetic rate constant (mol/L/s or similar units)
                
                !> ------------------------------------------------------------
                !> Chapuza: Special handling for aerobic DOC degradation
                !> ------------------------------------------------------------
            !> Chapuza
                if (this%redox_kin_reacts(react_ind)%name=='aerobic degradation of DOC (review)') then
                    !> Modify stoichiometry with yield coefficients for aerobic degradation
                    !> We multiply stoichiometry by yield coefficients
                    !! We assume that the order of species is the same as in review paper Carrera et al (2022), page 26
                    
                    !> Apply yield coefficient to species 2, 3, 4 (biomass production pathway)
                    this%redox_kin_reacts(react_ind)%stoichiometry(2:4)=this%redox_kin_reacts(react_ind)%stoichiometry(2:4)*&
                        this%redox_kin_reacts(react_ind)%yield
                        !> multiply coefficients by yield (fraction converted to biomass)
                    
                    !> Apply yield coefficient to species 6
                    this%redox_kin_reacts(react_ind)%stoichiometry(6)=this%redox_kin_reacts(react_ind)%stoichiometry(6)*&
                        this%redox_kin_reacts(react_ind)%yield
                        !> multiply coefficient by yield
                    
                    !> Apply (1-yield) to species 5 (energy/respiration pathway)
                    this%redox_kin_reacts(react_ind)%stoichiometry(5)=this%redox_kin_reacts(react_ind)%stoichiometry(5)*&
                        (1d0-this%redox_kin_reacts(react_ind)%yield)
                        !> multiply coefficient by (1-yield) (fraction used for energy)
                    
                    !> Adjust species 7 coefficient with yield correction
                    this%redox_kin_reacts(react_ind)%stoichiometry(7)=this%redox_kin_reacts(react_ind)%stoichiometry(7)-&
                        2d-1*this%redox_kin_reacts(react_ind)%yield
                        !> subtract 0.2*yield from coefficient (mass balance adjustment)
                end if
                !> Map species names to chemical system indices
                call this%redox_kin_reacts(react_ind)%set_species_indices_from_names(tmp_sp_names,this%species)
                
            else !> not a kinetic reaction - check if equilibrium reaction
                !> ----------------------------------------------------------------
                !> Process equilibrium reaction
                !> ----------------------------------------------------------------
                call this%is_eq_reaction_in_chem_syst(label,flag,react_ind)
                    !> Search for reaction in equilibrium reactions list
                    !> Returns: flag (TRUE if found), react_ind (index in system)
                
                if (flag.eqv. .true.) then !> if reaction is an equilibrium reaction
                    !> Read equilibrium reaction structural info
                    read(unit,*) num_inh !> read number of inhibitors (informational only)
                    read(unit,*) num_TEAs !> read number of terminal electron acceptors (informational only)
                    
                    !> Read reaction stoichiometry, species, rate constant, and equilibrium constant
                    if (allocated(tmp_sp_names)) deallocate(tmp_sp_names)
                    allocate(tmp_sp_names(this%eq_reacts(react_ind)%num_species))
                    read(unit,*) this%eq_reacts(react_ind)%num_species, (tmp_sp_names(j), &
                        this%eq_reacts(react_ind)%stoichiometry(j), j=1,this%eq_reacts(react_ind)%num_species), rate_cst, logK
                        !> reads: n_sp, name_1, coeff_1, ..., name_n, coeff_n, rate_cst, log10(K_eq)
                    
                    !> Set equilibrium constant (convert from log K)
                    call this%eq_reacts(react_ind)%set_eq_cst(10**logK)
                        !> K_eq = 10^(log K)
                else !> reaction not found in either kinetic or equilibrium list
                    continue !> skip this reaction, continue to next
                end if
            end if
        end if
    end do !> end second iteration loop
    
!> ================================================================
!> SECTION 5: Cleanup
!> ================================================================
!> Post-process
    close(unit) !> close Monod reactions database file
end subroutine