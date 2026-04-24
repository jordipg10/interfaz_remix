!> @file read_kinetics_DB.f90
!> @brief Reads kinetic mineral dissolution/precipitation parameters from database file
!>
!> @details
!> This subroutine reads kinetic mineral reaction parameters from the 'kinetics_modif.dat' database file.
!> The database contains information about kinetic mineral reactions including:
!> - Rate constants for parallel reactions
!> - Catalytic effects (proton-promoted, hydroxyl-promoted, etc.)
!> - Activation energies for temperature dependence
!> - Supersaturation thresholds for precipitation inhibition
!> - Stoichiometric coefficients for catalyzer species
!>
!> The kinetic rate law implemented follows the form:
!> \f[
!> r = \sum_{k=1}^{N_k} k_k \cdot \prod_{c=1}^{N_c} (a_c)^{p_{k,c}} \cdot (1 - \Omega^\theta)^\eta
!> \f]
!> where:
!> - \f$N_k\f$ = number of parallel reactions
!> - \f$k_k\f$ = rate constant for reaction k
!> - \f$N_c\f$ = number of catalyzers
!> - \f$a_c\f$ = activity of catalyzer c
!> - \f$p_{k,c}\f$ = power of catalyzer c in reaction k
!> - \f$\Omega\f$ = saturation index
!> - \f$\theta\f$ = Temkin coefficient
!> - \f$\eta\f$ = reaction order
!>
!> File format (each mineral):
!> ```
!> mineral_name  N_k  N_c  act_E  supersat_thr  comments
!> [For each parallel reaction k=1..N_k:]
!>   k_k  theta_k  eta_k
!>   [For each catalyzer c=1..N_c:]
!>     catalyzer_name  p_k_c
!> ```
!>
!> @param[in,out] this Chemical system object containing minerals and aqueous species
!> @param[in] path Directory path where 'kinetics_modif.dat' is located
!> @param[in] unit File unit number for reading
!>
!> @note The activation energy is read in kcal/mol and converted to J/mol (factor: 1484)
!> @note Only minerals with index <= num_minerals_kin are processed (kinetic minerals)
!> @note The file must end with a line containing 'null' as the mineral name
!>
!> @warning If a catalyzer species is not found in the aqueous phase, the program stops with error
!> @warning Minerals not in the chemical system are silently skipped
!>
!> @see kin_mineral_params_c
!> @see kin_mineral_c
!>
!> @author jordi Prats
!> @date 2025
!>
subroutine read_kinetics_DB(this,path,unit)
    !> ================================================================
    !> Module imports - bring in required chemistry types
    !> ================================================================
    use chem_system_m, only: chem_system_c
    use mineral_m, only: mineral_c
    use kin_mineral_m, only: kin_mineral_c
    use kin_mineral_params_m, only: kin_mineral_params_c
    use aq_species_m, only: aq_species_c

    implicit none !> enforce explicit variable declarations
    
    !> ================================================================
    !> Input/output parameters
    !> ================================================================
    class(chem_system_c) :: this                    !> chemical system object (modified - kinetic params added)
    character(len=*), intent(in) :: path            !> directory path to database file
    integer(kind=4), intent(in) :: unit             !> file unit number for reading
    
    !> ================================================================
    !> Local variables - counters and indices
    !> ================================================================
    integer(kind=4) :: i,j,k,l,N_k,N_c,mineral_ind,cat_ind
        !> i: counter for kinetic minerals in chemical system
        !> j: counter for parallel reactions (1 to N_k)
        !> k: counter for catalyzers (1 to N_c)
        !> l: (not currently used)
        !> N_k: number of parallel reactions for current mineral
        !> N_c: number of catalyzers for current mineral
        !> mineral_ind: index of current mineral in chemical system
        !> cat_ind: index of catalyzer species in aqueous phase
    
    !> Allocatable integer arrays
    integer(kind=4), allocatable :: global_stoich_counter(:),global_stoich_indices(:),aux_ind(:)
        !> global_stoich_counter: (not currently used) - for global stoichiometry
        !> global_stoich_indices: (not currently used) - species indices
        !> aux_ind: (not currently used) - auxiliary index array
    
    !> Logical flag
    logical :: flag !> TRUE if mineral/species found in chemical system, FALSE otherwise
    
    !> ================================================================
    !> Local variables - kinetic parameters
    !> ================================================================
    real(kind=8) :: unk,logK,valence,act_E,supersat_thr,cst
        !> unk: (not currently used)
        !> logK: (not currently used) - equilibrium constant
        !> valence: (not currently used) - species valence
        !> act_E: activation energy (read in kcal/mol, converted to J/mol)
        !> supersat_thr: supersaturation threshold for precipitation inhibition
        !> cst: constant (used when skipping non-kinetic minerals)
    
    !> Allocatable real arrays
    real(kind=8), allocatable :: global_stoich_coeffs(:),aux_coeffs(:)
        !> global_stoich_coeffs: (not currently used) - stoichiometric coefficients
        !> aux_coeffs: (not currently used) - auxiliary coefficients
    
    !> ================================================================
    !> Character variables for parsing
    !> ================================================================
    character(len=256) :: name,comments,label,filename
        !> name: mineral name read from file
        !> comments: comment field for mineral (not used)
        !> label: section label (first line of file)
        !> filename: full path to kinetics database file
    character(len=256), allocatable :: str1(:),str2(:),species_str(:),global_stoich_names(:)
        !> str1: (not currently used) - string array 1
        !> str2: (not currently used) - string array 2
        !> species_str: (not currently used) - species names
        !> global_stoich_names: (not currently used) - stoichiometry species names
    
    !> ================================================================
    !> Temporary chemistry objects
    !> ================================================================
    type(mineral_c) :: mineral !> temporary mineral object for searching
    type(kin_mineral_c), allocatable :: min_kin_reacts(:) !> array of kinetic mineral reactions (not used)
    type(kin_mineral_params_c) :: min_params !> kinetic parameters object for current mineral
    
    type(aq_species_c), allocatable :: catalysers(:) !> array of catalyzer species for current mineral
    
    !> ================================================================
    !> SECTION 1: Initialization and file opening
    !> ================================================================
    !> Allocate kinetic mineral reactions array (not currently used in this subroutine)
    allocate(min_kin_reacts(this%num_minerals_kin)) !> size = number of kinetic minerals in system
    
    !> Construct full filename path
    filename=trim(path)//'\kinetics_modif.dat' !> concatenate path + '\kinetics_modif.dat'

    !> Open kinetics database file for reading
    open(unit,file=filename,status='old',action='read') !> status='old' means file must exist
    
    !> Initialize counters
    i=1 !> counter for kinetic minerals in chemical system (incremented after each successful read)
    j=1 !> counter for parallel reactions (reset for each mineral)
    k=1 !> counter for catalysts (reset for each reaction)
    
    !> Read header label (first line of file)
    read(unit,*) label !> section label/header (not used)
    
    !> ================================================================
    !> SECTION 2: Main reading loop - process each mineral
    !> ================================================================
    do !> infinite loop - exit when 'null' encountered or all kinetic minerals processed
        !> ----------------------------------------------------------------
        !> Read mineral header line
        !> ----------------------------------------------------------------
        read(unit,*) mineral%name, N_k, N_c, act_E, supersat_thr, comments
            !> mineral%name: name of mineral
            !> N_k: number of parallel reactions for this mineral
            !> N_c: number of catalyzers (applies to all reactions)
            !> act_E: activation energy (kcal/mol)
            !> supersat_thr: supersaturation threshold
            !> comments: comment string (not used)
        
        !> ----------------------------------------------------------------
        !> Check exit condition
        !> ----------------------------------------------------------------
        if (name=='null') then !> if sentinel value found (end of data)
            exit !> exit reading loop
        else !> continue processing mineral
            !> ----------------------------------------------------------------
            !> Check if mineral exists in chemical system
            !> ----------------------------------------------------------------
            call this%is_mineral_in_chem_syst(mineral,flag,mineral_ind)
                !> Search for mineral in chemical system
                !> Returns: flag (TRUE if found), mineral_ind (index in system)
            
            !> ----------------------------------------------------------------
            !> Process mineral if it's a kinetic mineral
            !> ----------------------------------------------------------------
            if (flag.eqv..true. .and. mineral_ind<=this%num_minerals_kin) then
                !> Mineral found AND it's kinetic (index within kinetic range)
                
                !> ------------------------------------------------------------
                !> Set basic kinetic parameters
                !> ------------------------------------------------------------
                min_params%num_par_reacts=N_k !> store number of parallel reactions
                min_params%num_cat=N_c !> store number of catalyzers
                min_params%act_energy=act_E*1484 !> convert kcal/mol to J/mol (factor: 1484 = 4184/2.82)
                min_params%supersat_threshold=supersat_thr !> store supersaturation threshold
                
                !> Allocate parameter arrays
                call min_params%allocate_constants() !> allocate k, theta, eta arrays (size N_k)
                call min_params%allocate_cat_indices() !> allocate p array (N_k x N_c) and cat_indices (N_c)
                allocate(catalysers(min_params%num_cat)) !> allocate catalyzer species array (size N_c)
                
                !> ------------------------------------------------------------
                !> Read parallel reactions and their catalyzers
                !> ------------------------------------------------------------
                do j=1,min_params%num_par_reacts !> loop over each parallel reaction
                    !> Read reaction rate parameters
                    read(unit,*) min_params%k(j), min_params%theta(j), min_params%eta(j)
                        !> k(j): rate constant for reaction j
                        !> theta(j): Temkin coefficient for reaction j
                        !> eta(j): reaction order for reaction j
                    
                    !> Read catalyzer information for this reaction
                    do k=1,min_params%num_cat !> loop over each catalyzer
                        read(unit,*) catalysers(k)%name, min_params%p(j,k)
                            !> catalysers(k)%name: name of catalyzer species k
                            !> p(j,k): power of catalyzer k in reaction j
                    end do
                    
                    !> Find catalyzer indices in chemical system
                    do k=1,min_params%num_cat !> loop over each catalyzer again
                        !> Search for catalyzer in aqueous phase
                        call this%aq_phase%is_species_in_aq_phase(catalysers(k),flag,cat_ind)
                            !> Returns: flag (TRUE if found), cat_ind (index in aqueous phase)
                        
                        if (flag.eqv..false.) then !> if catalyzer not found
                            error stop 'Catalyser not found in chemical system' !> stop with error
                        else !> catalyzer found
                            min_params%cat_indices(k)=cat_ind !> store catalyzer index
                        end if
                    end do
                end do !> end parallel reactions loop
                
                !> ------------------------------------------------------------
                !> Store parameters in chemical system
                !> ------------------------------------------------------------
                call this%min_kin_reacts(mineral_ind)%set_mineral_params(min_params)
                    !> Assign kinetic parameters to mineral at index mineral_ind
                call this%min_kin_reacts(mineral_ind)%params%set_m(9.3d-1) !> set m parameter for kinetic rate law (can be adjusted as needed)
                !> ------------------------------------------------------------
                !> Check if more minerals to process
                !> ------------------------------------------------------------
                if (i<this%num_minerals_kin) then !> if not all kinetic minerals processed yet
                    i=i+1 !> increment kinetic mineral counter
                else !> all kinetic minerals processed
                    exit !> exit reading loop
                end if
                
                !> Deallocate temporary arrays for next mineral
                deallocate(min_params%k,min_params%p,min_params%theta,min_params%eta,min_params%cat_indices,catalysers)
                    !> Free memory before next iteration
                
            else !> mineral not found OR not kinetic
                !> ----------------------------------------------------------------
                !> Skip mineral data (read and discard)
                !> ----------------------------------------------------------------
                read(unit,*) cst !> read one dummy line to skip mineral data
            end if
        end if
    end do !> end main reading loop
    
    !> ================================================================
    !> SECTION 3: Cleanup and file closing
    !> ================================================================
    close(unit) !> close kinetics database file
    
end subroutine !> end of read_kinetics_DB