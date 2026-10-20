!> \file write_chemistry.f90
!> \brief Writes comprehensive chemistry simulation results to formatted output file
!> \details This subroutine writes all chemical system configuration and simulation results
!> to a formatted text file, including:
!> - Chemical system composition (species, minerals, gases, reactions)
!> - Stoichiometric matrices
!> - Concentrations and activities (aqueous, solid, gas phases)
!> - Reaction rates and amounts (equilibrium and kinetic)
!> - Zone-specific data (reactive zones, mineral zones, cation exchange zones)
!> - Lagrangian particle information (if applicable)
!> - Water types (recharge, boundary, domain)
!> - pH, salinity, and ionic strength
!> - Volumetric fractions and partial pressures
!> 
!> Output format uses scientific notation (ES15.5) for real values and fixed-width formatting
!> for species names. Data is organized by zone type with clear section headers.

!> \brief Write chemistry object data and results to output file
!> \details Comprehensive output subroutine that writes all chemical system information including:
!> chemical system configuration, species lists, reactions, concentrations, activities, reaction rates,
!> and zone-specific data. Handles both Eulerian and Lagrangian formulations.
!> 
!> \param[in] this chemistry_c object containing all chemical system data and results
!> \param[in] unit Integer file unit number for output [-]
!> \param[in] Lagr_flag Logical flag indicating Lagrangian formulation (true) or Eulerian (false) [-]
!> 
!> \use chemistry_m Chemistry module providing chemistry_c class definition
!> Writes data and results of chemistry object
subroutine write_chemistry(this,dir,root,time_step,Lagr_flag,num_time)
    use chemistry_m, only: chemistry_c  !> Import chemistry class and infinity norm function
    use vectors_m, only: inf_norm_vec_real !> Import vector class for real numbers
    implicit none
    class(chemistry_c), intent(in) :: this !> chemistry object containing all simulation data
    character(len=*), intent(in) :: dir !> directory for output files
    character(len=*), intent(in) :: root !> root name for output files
    integer(kind=4), intent(in), optional :: time_step !> time step index for output file naming [-]
    logical, intent(in), optional :: Lagr_flag !> transport approach: .true.=Lagrangian, .false.=Eulerian (used in output filename) [-]
    integer(kind=4), intent(in), optional :: num_time !> total number of time steps in the run (appended to output filename) [-]
    
    integer(kind=4) :: i,j,l,num_wat,num_ext !> loop indices and counters [-]
    integer(kind=4), allocatable :: wat_indices(:),tar_sol_indices(:) !> domain and external water indices [-]
    real(kind=8), allocatable :: sat_row(:) !> saturation index of one kinetic mineral in each target water (1 = equilibrium where never evaluated) [-]
    integer(kind=4) :: unit !> file unit number for main output [-]
    integer(kind=4) :: jw !> water index for species order check [-]
    logical :: order_ok !> flag for species order consistency check [-]
    character(len=:), allocatable :: wma_tag !> WMA option tag for output filename (lumped/consistent) [-]
    character(len=:), allocatable :: tpt_tag !> transport approach tag for output filename (Euler/Lagr) [-]
    character(len=:), allocatable :: nt_tag !> number-of-time-steps tag for output filename [-]
    character(len=32) :: nt_buf !> scratch buffer for formatting the time-step count [-]
    logical :: mesh_too_large !> .true. if the mesh is too large to dump per-cell data to the text file [-]
    integer(kind=4), parameter :: max_mesh_print = 200 !> per-cell concentrations/reaction rates are written to the text file only up to this many target waters; larger meshes go to the binary file [-]
    
    if (this%lump_flag .eqv. .true.) then !> select WMA option tag for filename
        wma_tag = '_lumped' !> lumped WMA
    else
        wma_tag = '_consistent' !> consistent WMA
    end if
    if (present(Lagr_flag)) then !> select transport approach tag for filename
        if (Lagr_flag .eqv. .true.) then
            tpt_tag = '_Lagr' !> Lagrangian transport
        else
            tpt_tag = '_Euler' !> Eulerian transport
        end if
    else
        tpt_tag = '' !> transport approach unspecified: omit tag
    end if
    if (present(num_time)) then !> append the total number of time steps to the filename
        write(nt_buf,'(I0)') num_time
        nt_tag = '_'//trim(nt_buf)//'steps'
    else
        nt_tag = '' !> number of time steps unspecified: omit tag
    end if
    open(newunit=unit, file=dir//root//tpt_tag//wma_tag//nt_tag//'_chemistry.out', status='unknown', form='formatted')
    !> For large meshes, the per-cell concentration/activity/reaction blocks bloat the text file;
    !> they are skipped here (a compact binary file written by write_mesh_bin holds them instead).
    mesh_too_large = this%num_target_waters > max_mesh_print
    write(unit,"(2x,'Chemical input option:',I25)") this%read_opt !> write chemical input option flag
    if (this%lump_flag .eqv. .true. ) then !> check if lumped WMA is used
        write(unit,"(/,2x,'WMA option:',10x,'Lumped')") !> write lumped WMA label
    else !> consistent WMA is used
        write(unit,"(/,2x,'WMA option:',10x,'Consistent')") !> write consistent WMA label
        !write(unit,"(/,2x,'Consistent WMA option:',I25)") this%cons_opt !> write consistent WMA option flag
        write(unit,"(/,2x,'rk downstream waters option:',I25)") this%r_down_opt !> write kinetic reaction downstream option
        if (this%r_down_opt>1) then !> check if coefficient of variation parameter is used
            write(unit,"(/,2x,'Parameter mu:',ES15.5)") this%CV_params%est_prm !> write CV parameter mu [-]
        end if
    end if
    write(unit,"(/,2x,'rk average option:',I25)") this%rk_avg_opt !> write kinetic reaction averaging option flag
    write(unit,"(/,2x,'Aqueous species:'/)")  !> write aqueous species section header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,A15)") this%chem_syst%aq_phase%aq_species(i)%name !> write aqueous species name
    end do
    write(unit,"(2x,'Minerals (name + molar volume):'/)")  !> write minerals section header
    do i=1,this%chem_syst%num_minerals !> loop over all minerals in chemical system
        write(unit,"(10x,A15,ES15.5)") this%chem_syst%minerals(i)%name, this%chem_syst%minerals(i)%mineral%mol_vol !> write mineral name and molar volume [dm³/mol]
    end do
    write(unit,"(2x,'Gases:'/)")  !> write gases section header
    do i=1,this%chem_syst%gas_phase%num_species !> loop over all gas species
        write(unit,"(10x,A15)") this%chem_syst%gas_phase%gases(i)%name !> write gas species name
    end do
    write(unit,"(2x,'Surface complexes:'/)")  !> write surface complexes section header
    do i=1,this%chem_syst%cat_exch_zone%num_surf_compl !> loop over all surface complexation species
        write(unit,"(10x,A15)") this%chem_syst%cat_exch_zone%surf_compl(i)%name !> write surface complex name
    end do
    write(unit,"(2x,'Equilibrium reactions:'/)")  !> write equilibrium reactions section header
    do i=1,this%chem_syst%speciation_alg%num_eq_reactions !> loop over all equilibrium reactions
        write(unit,"(10x,A30/)") this%chem_syst%eq_reacts(i)%name !> write equilibrium reaction name
    end do
    write(unit,"(2x,'Kinetic reactions:'/)")  !> write kinetic reactions section header
    do i=1,this%chem_syst%num_lin_kin_reacts
        call this%chem_syst%lin_kin_reacts(i)%write_reaction(unit)
    end do
    do i=1,this%chem_syst%num_redox_kin_reacts
        call this%chem_syst%redox_kin_reacts(i)%write_reaction(unit)
    end do
    do i=1,this%chem_syst%num_minerals_kin
        call this%chem_syst%min_kin_reacts(i)%write_reaction(unit)
    end do
    write(unit,"(2x,'Global stoichiometric matrix:',/)")  !> write stoichiometric matrix section header
    do i=1,this%chem_syst%num_reacts !> loop over all reactions (equilibrium + kinetic)
        write(unit,"(10x,*(F15.5))") (this%chem_syst%stoich_mat(i,j), j=1,this%chem_syst%speciation_alg%num_species) !> write stoichiometric coefficients for reaction i [-]
    end do
    write(unit,"(2x,'Equilibrium constants (database reference values, 25 C; see the water types block for their values at each temperature):',/)")  !> write equilibrium constants section header
    do i=1,this%chem_syst%speciation_alg%num_eq_reactions !> loop over all equilibrium reactions
        write(unit,"(10x,*(ES15.5))") this%chem_syst%eq_reacts(i)%eq_cst_ref !> write equilibrium constant K_eq [-]
    end do
    write(unit,"(/,2x,'Global component matrix:'/)")  !> write component matrix section header (currently commented out)
    ! do i=1,this%chem_syst%speciation_alg%num_prim_species
    !     write(unit,"(10x,*(F15.5))") (this%chem_syst%speciation_alg%comp_mat(i,j), &
    !     j=1,this%chem_syst%speciation_alg%num_var_act_species)
    ! end do
    !if (Lagr_flag) then !> check if Lagrangian formulation is used
        ! write(unit,"(2x,'Initial water particles coordinates (by rows):'/)")
        ! do i=1,this%num_waters
        !     write(unit,"(2x,*(ES15.5))") this%waters_init(i)%solid_chemistry%tar%coord !> Coordinates of the targets
        ! end do
        ! write(unit,"(/,2x,'Final water particles coordinates:'/)")
        ! do i=1,this%num_waters
        !     write(unit,"(2x,*(ES15.5))") this%waters(i)%solid_chemistry%tar%coord !> Coordinates of the targets
        ! end do
        !write(unit,"(2x,'Initial water particles positions:'/)")  !> write initial Lagrangian particle positions header
        !do i=1,this%num_target_waters_init !> loop over initial target waters
        !    write(unit,"(2x,I5,ES15.5)") this%waters_init(&
        !        this%tar_wat_indices_init(i))%id, &
        !        this%waters_init(this%tar_wat_indices_init(i))%pos(1) !> write particle ID and spatial position [L]
        !end do
        !write(unit,"(/,2x,'Final water particles positions:'/)")  !> write final Lagrangian particle positions header
        !do i=1,this%num_target_waters !> loop over final target waters
        !    write(unit,"(2x,I5,ES15.5)") this%waters(this%tar_wat_indices(i))%id, &
        !        this%waters(this%tar_wat_indices(i))%pos(1) !> write particle ID and final position [L]
        !end do
        ! write(unit,"(/,2x,'Initial water particles volumes (by rows):'/)")  !> write initial particle volumes header
        ! do i=1,this%num_target_waters_init !> loop over initial target waters
        !     write(unit,"(2x,ES15.5)") this%waters_init(i)%volume !> write initial particle volume [L³]
        ! end do
        ! write(unit,"(/,2x,'Final water particles volumes:'/)")  !> write final particle volumes header
        ! do i=1,this%num_target_waters !> loop over final target waters
        !     write(unit,"(2x,ES15.5)") this%waters(i)%volume !> write final particle volume [L³]
        ! end do
    !end if
    ! write(unit,"(/,2x,'Concentration of recharge waters (in molalities): &
    !     (rows -> aqueous species in the same order as before, columns -> targets)'/)") 
    !     !> write recharge water concentrations header
    ! do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
    !     write(unit,"(10x,*(ES15.5))") (this%waters(this%rech_waters_indices(j))%concentrations(&
    !         this%waters(this%rech_waters_indices(j))%indices_aq_species(i)), j=1,this%num_rech_waters) !> write recharge water concentrations [C]
    ! end do
    ! write(unit,"(/,2x,'Activities of recharge waters:'/)")  !> write recharge water activities header
    ! do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
    !     write(unit,"(10x,*(ES15.5))") (this%waters(this%rech_waters_indices(j))%activities(&
    !         this%waters(this%rech_waters_indices(j))%indices_aq_species(i)), j=1,this%num_rech_waters) !> write recharge water activities [-]
    ! end do
    !> Check that all water types have the same species ordering
    order_ok = .true.
    do j = 2, this%num_wat_types
        do i = 1, this%chem_syst%aq_phase%num_species
            if (this%wat_types(j)%indices_aq_species(i) /= this%wat_types(1)%indices_aq_species(i)) then
                print *, "WARNING: water type", j, " has different indices_aq_species than water type 1"
                print *, "  species", i, ": water type 1 ->", this%wat_types(1)%indices_aq_species(i), &
                    ", water type", j, " ->", this%wat_types(j)%indices_aq_species(i)
                order_ok = .false.
            end if
        end do
    end do
    !> Check that all target waters have the same species ordering as last water type
    do jw = 1, this%num_target_waters
        j = this%tar_wat_indices(jw)
        do i = 1, this%chem_syst%aq_phase%num_species
            if (this%waters(j)%indices_aq_species(i) /= this%wat_types(this%num_wat_types)%indices_aq_species(i)) then
                print *, "WARNING: target water", jw, " (index", j, ") has different indices_aq_species"
                print *, "  species", i, ": water type", this%num_wat_types, "->", &
                    this%wat_types(this%num_wat_types)%indices_aq_species(i), &
                    ", target water ->", this%waters(j)%indices_aq_species(i)
                order_ok = .false.
            end if
        end do
    end do
    !> Check initial target waters too
    do jw = 1, this%num_target_waters_init
        j = this%tar_wat_indices_init(jw)
        do i = 1, this%chem_syst%aq_phase%num_species
            if (this%waters_init(j)%indices_aq_species(i) /= this%wat_types(1)%indices_aq_species(i)) then
                print *, "WARNING: init target water", jw, " (index", j, ") has different indices_aq_species"
                print *, "  species", i, ": water type 1 ->", this%wat_types(1)%indices_aq_species(i), &
                    ", init target water ->", this%waters_init(j)%indices_aq_species(i)
                order_ok = .false.
            end if
        end do
    end do
    if (order_ok) then
        print *, "Species ordering check PASSED: all waters have consistent indices_aq_species"
    end if
    write(unit,"(/,2x,'Water types names:'/)")  !> write water types names header
    write(unit,"(/,10x,*(A25))") (this%wat_types(j)%name, j=1,this%num_wat_types)  !> write water types names (column width = ES25.15E3 below)
    !> Per water type: temperature, Debye-Huckel parameters and the equilibrium constants of its
    !> reactive zone AT ITS TEMPERATURE (the constants its speciation satisfies; the 'Equilibrium
    !> constants' blocks of the chemical system and of the reactive zones are the database
    !> reference values at 25 C). Each water type has its own reactive zone, hence its own set
    !> and order of equilibrium reactions: one sub-block per water type, reactions by name.
    write(unit,"(/,2x,'Equilibrium constants of each water type at its temperature (K(T), log10 K(T)) next to the database reference value (log10 K at 25 C):')")
    do j=1,this%num_wat_types !> loop over water types
        write(unit,"(/,10x,A,':   T = ',F8.3,' C   rho_w = ',F8.5,' kg/L   Debye-Huckel A = ',F8.5,'   B = ',F8.5)") &
            trim(this%wat_types(j)%name), this%wat_types(j)%temp-273.15d0, this%wat_types(j)%density, &
            this%wat_types(j)%params_aq_sol%A, this%wat_types(j)%params_aq_sol%B
        if (allocated(this%wat_types(j)%logK) .and. associated(this%wat_types(j)%solid_chemistry)) then
            write(unit,"(10x,A25,3A18)") 'reaction', 'K(T)', 'log10 K(T)', 'log10 K(25 C)'
            do i=1,this%wat_types(j)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !> equilibrium reactions of its reactive zone
                associate(rct => this%wat_types(j)%solid_chemistry%reactive_zone%chem_syst%eq_reacts(&
                    this%wat_types(j)%solid_chemistry%reactive_zone%ind_eq_reacts(i)))
                if (rct%eq_cst_ref>0d0) then
                    write(unit,"(10x,A25,ES18.5,2F18.6)") trim(rct%name), 1d1**this%wat_types(j)%logK(i), &
                        this%wat_types(j)%logK(i), log10(rct%eq_cst_ref)
                else
                    write(unit,"(10x,A25,ES18.5,F18.6)") trim(rct%name), 1d1**this%wat_types(j)%logK(i), &
                        this%wat_types(j)%logK(i)
                end if
                end associate
            end do
        end if
    end do
    write(unit,"(/,2x,'Concentration of water types (in molalities):'/)")  !> write water types concentrations header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES25.15E3))") (this%wat_types(j)%concentrations(&
            this%wat_types(j)%indices_aq_species(i)), j=1,this%num_wat_types) !> write water type concentrations [C] (full double precision: 16 significant digits)
    end do
    write(unit,"(/,2x,'Activities of water types (in molalities):'/)")  !> write water types activities header
    !> Write water types concentrations header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES25.15E3))") (this%wat_types(j)%activities(&
            this%wat_types(j)%indices_aq_species(i)), j=1,this%num_wat_types) !> write water type activities [-] (full double precision: 16 significant digits)
    end do
    !if (this%num_bd_waters>0) then !> check if boundary waters are defined
        ! write(unit,"(/,2x,'Activities of boundary waters:'/)")  !> write boundary water activities header
        ! do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        !     write(unit,"(10x,*(ES15.5))") (this%waters(this%bd_waters_indices(j))%activities(&
        !         this%waters(this%bd_waters_indices(j))%indices_aq_species(i)), j=1,this%num_bd_waters) !> write boundary water activities [-]
        ! end do
    !end if
    if (mesh_too_large) then
        write(unit,"(/,2x,'NOTE: mesh has',I10,' target waters (> ',I0,'): per-cell concentrations,')") &
            this%num_target_waters, max_mesh_print
        write(unit,"(2x,'activities and reaction rates/amounts are written to the binary mesh file')")
        write(unit,"(2x,'(',A,'_mesh.bin + _meta.txt) instead of this text file.'/)") trim(root)//tpt_tag//wma_tag
    else
    write(unit,"(/,2x,'Initial concentration of aqueous species in target waters:',/ &
        2x,'(rows -> aqueous species in the same order as before, columns -> target waters)'/)")  !> write initial domain concentrations header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices_init(j))%concentrations(&
            this%waters_init(this%tar_wat_indices_init(j))%indices_aq_species(i)), j=1,this%num_target_waters_init) !> write initial target water concentrations [C]
    end do
    write(unit,"(/,2x,'Initial activities aqueous species target waters:'/)")  !> write initial domain activities header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices_init(j))%activities(&
            this%waters_init(this%tar_wat_indices_init(j))%indices_aq_species(i)), j=1,this%num_target_waters_init) !> write initial target water activities [-]
    end do
    write(unit,"(/,2x,'Final concentration of aqueous species in target waters:'/)")  !> write final target water concentrations header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%concentrations(& 
            this%waters(this%tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters) !> write final target water concentrations [C]
    end do
    ! write(unit,"(/,2x,'Final concentration of aqueous species in the mesh:'/)")  !> write final mesh concentrations header
    ! do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
    !    write(unit,"(10x,*(ES15.5))") this%waters(this%bd_waters_indices(1))%concentrations(this%waters(&
    !        this%bd_waters_indices(1))%indices_aq_species(i)), &
    !        (this%waters(this%tar_wat_indices(j))%concentrations(& 
    !        this%waters(this%tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters) & 
    !        , this%waters(this%bd_waters_indices(this%num_bd_waters))%concentrations(this%waters( &
    !        this%bd_waters_indices(this%num_bd_waters))%indices_aq_species(i)) !> write final mesh concentrations [C]
    ! end do
    write(unit,"(/,2x,'Final activities aqueous species target waters:'/)")  !> write final domain activities header
    do i=1,this%chem_syst%aq_phase%num_species !> loop over all aqueous species
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%activities(&
            this%waters(this%tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters) !> write final domain activities [-]
    end do
    !if (this%chem_syst%aq_phase%ind_wat>0) then
        ! write(unit,"(/,2x,'Initial activity water domain:'/)")
        ! write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices(j))%activities(&
        !     this%waters_init(this%tar_wat_indices(j))%indices_aq_species(this%chem_syst%aq_phase%ind_wat)),&
        !     j=1,this%num_target_waters)
        write(unit,"(/,2x,'Initial salinity target waters:'/)")  !> write initial domain salinity header
        write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices(j))%salinity,j=1,this%num_target_waters) !> write initial salinity [C]
        ! write(unit,"(/,2x,'Activity water domain:'/)")
        ! write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%activities(&
        !     this%waters(this%tar_wat_indices(j))%indices_aq_species(this%chem_syst%aq_phase%ind_wat)),&
        !     j=1,this%num_target_waters)
        write(unit,"(/,2x,'Final salinity target waters:'/)")  !> write final domain salinity header
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%salinity, &
            j=1,this%num_target_waters) !> write final salinity [C]
    !end if
    if (this%chem_syst%aq_phase%ind_prot>0) then !> check if proton is present in system
        write(unit,"(/,2x,'Initial pH target waters:'/)")  !> write initial pH header
        write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices(j))%pH, j=1,this%num_target_waters) !> write initial pH [-]
        write(unit,"(/,2x,'Final pH target waters:'/)")  !> write final pH header
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%pH, j=1,this%num_target_waters) !> write final pH [-]
    end if
    write(unit,"(/,2x,'Aqueous equilibrium reaction amounts:'/)")  !> write aqueous equilibrium reaction amounts header
    do i=1,this%chem_syst%aq_phase%num_aq_complexes !> loop over all aqueous complexes
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%rxn_amnt_eq(i), j=1,this%num_target_waters) !> write equilibrium reaction amount R_e [mol]
    end do
    write(unit,"(/,2x,'Aqueous equilibrium reaction rates:'/)")  !> write aqueous equilibrium reaction rates header
    do i=1,this%chem_syst%aq_phase%num_aq_complexes
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%rxn_rate_eq(i), j=1,this%num_target_waters)
    end do
    write(unit,"(/,2x,'Initial aqueous kinetic reaction rates:'/)")  !> write initial aqueous kinetic rates header
    do i=1,this%chem_syst%num_aq_kin_reacts !> loop over all aqueous kinetic reactions
        write(unit,"(10x,*(ES15.5))") (this%waters_init(this%tar_wat_indices(j))%rxn_amnt_kin(i), j=1,this%num_target_waters) !> write initial kinetic rate r_k [1/T]
    end do
    
    write(unit,"(/,2x,'Aqueous mean kinetic reaction rates:'/)")  !> write mean aqueous kinetic rates header
    do i=1,this%chem_syst%num_aq_kin_reacts !> loop over all aqueous kinetic reactions
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%rxn_rate_kin_avg(i), j=1,this%num_target_waters) !> write mean kinetic rate r_k_mean [1/T]
    end do
    write(unit,"(/,2x,'Aqueous kinetic reaction amounts:'/)")  !> write aqueous kinetic reaction amounts header
    do i=1,this%chem_syst%num_aq_kin_reacts !> loop over all aqueous kinetic reactions
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%rxn_amnt_kin(i), j=1,this%num_target_waters) !> write kinetic reaction amount R_k [mol]
    end do
    write(unit,"(/,2x,'Aqueous estimated kinetic reaction amounts:'/)")  !> write estimated aqueous kinetic amounts header
    do i=1,this%chem_syst%num_aq_kin_reacts !> loop over all aqueous kinetic reactions
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%rxn_amnt_kin(i), j=1,this%num_target_waters) !> write kinetic amount R_k [mol]
    end do
    write(unit,"(/,2x,'Aqueous accumulated kinetic reaction amounts:'/)")  !> write accumulated aqueous kinetic amounts header
    do i=1,this%chem_syst%num_aq_kin_reacts !> loop over all aqueous kinetic reactions
        write(unit,"(10x,*(ES15.5))") (this%waters(this%tar_wat_indices(j))%cum_rxn_amnt_kin(i), j=1,this%num_target_waters) !> write accumulated kinetic amount R_k_accum [mol]
    end do
    end if !> mesh_too_large: end of per-cell aqueous concentration/activity/reaction blocks
    do i=1,this%num_init_cat_exch_zones !> loop over all initial cation exchange zones
        write(unit,"(/,2x,'Initial cation exchange site',I5,':'/)") i !> write cation exchange zone index
        write(unit,"(/,10x,'Initial concentration of surface complexes:'/)")  !> write initial surface complex concentrations header
        do j=1,this%init_cat_exch_zones(i)%reactive_zone%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
            write(unit,"(10x,*(ES25.15))") this%init_cat_exch_zones(i)%concentrations(j) !> write initial surface complex concentration [C]
        end do
        write(unit,"(/,10x,'Initial activity of surface complexes:'/)")  !> write initial surface complex activities header
        do j=1,this%init_cat_exch_zones(i)%reactive_zone%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
            write(unit,"(10x,*(ES25.15))") this%init_cat_exch_zones(i)%activities(j) !> write initial surface complex activity [-]
        end do
    end do
    if (.not. mesh_too_large .and. this%chem_syst%num_minerals_kin>0) then !> per-cell kinetic mineral saturation indices (skipped for large meshes like the other per-cell blocks)
        write(unit,"(/,2x,'Kinetic mineral saturation indices (IAP/K):'/)")  !> write kinetic mineral saturation header
        allocate(sat_row(this%num_target_waters))
        do i=1,this%chem_syst%num_minerals_kin !> loop over the kinetic minerals of the chemical system (they come first in chem_syst%minerals; the equilibrium-mineral entries of the saturation array are not written)
            do j=1,this%num_target_waters
                sat_row(j)=1d0 !> default: saturation never evaluated for this water/mineral (equilibrium, the attribute's initial value)
                if (allocated(this%waters(this%tar_wat_indices(j))%saturation_mins)) then
                    if (i<=size(this%waters(this%tar_wat_indices(j))%saturation_mins)) then
                        sat_row(j)=this%waters(this%tar_wat_indices(j))%saturation_mins(i)
                    end if
                end if
            end do
            write(unit,"(10x,A15,*(ES15.5))") this%chem_syst%minerals(i)%name, sat_row !> write mineral name and saturation index Omega per target water [-]
        end do
        deallocate(sat_row)
    end if
    do l=1,this%num_mineral_zones !> loop over all mineral zones
        write(unit,"(/,2x,'Mineral zone',I5,':'/)") l !> write mineral zone index
        write(unit,"(10x,'Minerals in equilibrium:'/)")  !> write equilibrium minerals section header
        do i=1,this%mineral_zones(l)%num_minerals_eq !> loop over equilibrium minerals in zone
            write(unit,"(20x,A15/)") this%mineral_zones(l)%chem_syst%minerals(&
                this%mineral_zones(l)%ind_min_chem_syst(this%mineral_zones(l)%num_minerals_kin+i))%name !> write equilibrium mineral name
        end do
        call this%link_target_solids_mineral_zone(l,tar_sol_indices) !> link target solids to current mineral zone
        !write(unit,"(/,10x,'Mineral equilibrium reaction rates:'/)")
        !do i=1,this%mineral_zones(l)%num_minerals_eq
        !    write(unit,"(10x,*(ES15.5))") (this%waters(wat_indices(j))%solid_chemistry%rxn_rate_eq(i), j=1,size(wat_indices))
        !end do
        write(unit,"(10x,'Kinetic minerals:'/)")  !> write kinetic minerals section header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals in zone
            write(unit,"(20x,A15/)")this%mineral_zones(l)%chem_syst%minerals(this%mineral_zones(l)%ind_min_chem_syst(i))%name !> write kinetic mineral name
        end do
        if (.not. mesh_too_large) then !> skip per-cell mineral rates/amounts for large meshes (see binary mesh file)
        write(unit,"(/,10x,'Mineral initial kinetic reaction rates:'/)")  !> write mineral kinetic rates header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            !print *, this%waters(wat_indices(1))%solid_chemistry%Rk_mean
            write(unit,"(10x,*(ES15.5))") (this%target_solids_init(tar_sol_indices(j))%rxn_rate_kin(i), j=1,SIZE(tar_sol_indices)) !> write mean mineral kinetic rate r_k_mean [1/T]
        end do
        write(unit,"(/,10x,'Mineral final kinetic reaction rates:'/)")  !> write mineral kinetic rates header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            !print *, this%waters(wat_indices(1))%solid_chemistry%Rk_mean
            write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_rate_kin(i), j=1,SIZE(tar_sol_indices)) !> write mean mineral kinetic rate r_k_mean [1/T]

        end do
        write(unit,"(/,10x,'Mineral mean kinetic reaction rates:'/)")  !> write mineral kinetic rates header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            !print *, this%waters(wat_indices(1))%solid_chemistry%Rk_mean
            write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_rate_kin_avg(i), j=1,SIZE(tar_sol_indices)) !> write mean mineral kinetic rate r_k_mean [1/T]
        end do
        write(unit,"(/,10x,'Mineral kinetic reaction amounts:'/)")  !> write mineral kinetic amounts header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            !print *, this%waters(wat_indices(1))%solid_chemistry%Rk_mean
            write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_amnt_kin(i), j=1,SIZE(tar_sol_indices)) !> write mineral kinetic amount R_k [mol]
        end do
        write(unit,"(/,10x,'Mineral estimated kinetic reaction amounts:'/)")  !> write estimated mineral kinetic amounts header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_amnt_kin(i), j=1,SIZE(tar_sol_indices)) !> write mineral kinetic amount R_k [mol]
        end do
        
        write(unit,"(/,10x,'Accumulated mineral kinetic reaction amounts [kg] (amount [mol] times the molecular weight of the mineral; 0 if the database gives none):'/)")  !> write accumulated mineral kinetic amounts header
        do i=1,this%mineral_zones(l)%num_minerals_kin !> loop over kinetic minerals
            write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%cum_rxn_amnt_kin(i)*&
                this%mineral_zones(l)%chem_syst%minerals(this%mineral_zones(l)%ind_min_chem_syst(i))%mineral%molecular_weight, &
                j=1,SIZE(tar_sol_indices)) !> write accumulated mineral kinetic amount R_k_accum [mol] x MW [kg/mol]
        end do
        end if !> mesh_too_large: end of per-cell mineral kinetic rates/amounts
    end do
    do l=1,this%num_reactive_zones !> loop over all reactive zones
        call this%link_target_solids_reactive_zone(l,tar_sol_indices) !> link target solids to current reactive zone
        call this%link_waters_target_solids(tar_sol_indices,wat_indices) !> link target waters to current reactive zone
        num_wat=size(wat_indices) !> get number of target waters in zone
        !num_ext=size(ext_indices) !> get number of external waters in zone
        if (num_wat>0) then !> check if zone has domain waters
            write(unit,"(/,2x,'Reactive zone',I5,':'/)") l !> write reactive zone index
            write(unit,"(10x,'Non-flowing species:'/)")  !> write non-flowing species section header
            do i=1,this%reactive_zones(l)%num_non_flow_species !> loop over non-flowing species in zone
                write(unit,"(20x,A15)") this%chem_syst%species(this%reactive_zones(l)%ind_non_flow_species(i))%name !> write non-flowing species name
            end do
            write(unit,"(10x,'Minerals:'/)")  !> write minerals section header
            do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in zone
                write(unit,"(20x,A15)") this%chem_syst%minerals(this%reactive_zones(l)%ind_mins_chem_syst(i))%name !> write mineral name
            end do
            write(unit,"(10x,'Number of equilibrium reactions:',I10/)") this%reactive_zones(l)%speciation_alg%num_eq_reactions !> write number of equilibrium reactions in zone
            write(unit,"(10x,'Se:',/)")  !> write equilibrium stoichiometric matrix header
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions !> loop over equilibrium reactions
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%stoich_mat(i,j), &
                j=1,this%reactive_zones(l)%speciation_alg%num_species) !> write equilibrium stoichiometric coefficients S_e [-]
            end do
            write(unit,"(10x,'Sk:',/)")  !> write kinetic stoichiometric matrix header
            do i=1,this%reactive_zones(l)%chem_syst%num_kin_reacts !> loop over kinetic reactions
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%chem_syst%Sk(i,j), &
                j=1,this%reactive_zones(l)%chem_syst%speciation_alg%num_species) !> write kinetic stoichiometric coefficients S_k [-]
            end do
            write(unit,"(10x,'U*Sk^T:',/)")  !> write U*Sk^T product matrix header
            if (allocated(this%reactive_zones(l)%U_SkT_prod)) then
                do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species !> loop over aqueous primary species
                    write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%U_SkT_prod(i,j), j=1,size(this%reactive_zones(l)%U_SkT_prod,2)) !> write U*Sk^T product matrix elements [-]
                end do
            end if
            write(unit,"(/,10x,'Equilibrium constants (database reference values, 25 C):',/)")  !> write zone equilibrium constants header
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions !> loop over equilibrium reactions in zone
                write(unit,"(10x,*(ES15.5))") this%reactive_zones(l)%chem_syst%eq_reacts(&
                    this%reactive_zones(l)%ind_eq_reacts(i))%eq_cst_ref !> write equilibrium constant K_eq [-]
            end do
            !write(unit,"(/,10x,'Aqueous species:',/)")
            !do i=1,this%chem_syst%aq_phase%num_species
            !    write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%aq_phase%aq_species(this%waters(wat_indices(1))%indices_aq_species(i))%name
            !end do
            !write(unit,"(/,10x,'Primary species:',/)")
            !do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species
            !    ! write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%species(&
            !    !     this%reactive_zones(l)%ind_aq_var_act_species(i))%name
            !end do
            if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl>0) then !> check if cation exchange zone has surface complexes
                write(unit,"(20x,A15)") this%reactive_zones(l)%cat_exch_zone%surf_compl(1)%name !> write first surface complex name
            end if
            write(unit,"(/,10x,'Component matrix:'/)")  !> write component matrix section header
            !print *, l
            !print *, this%reactive_zones(l)%speciation_alg%num_prim_species
            !print *, this%reactive_zones(l)%speciation_alg%num_var_act_species
            !print *, this%reactive_zones(l)%speciation_alg%comp_mat
            if (allocated(this%reactive_zones(l)%speciation_alg%comp_mat)) then !> check if component matrix is allocated
                do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species !> loop over primary species
                    write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%speciation_alg%comp_mat(i,j), &
                        j=1,this%reactive_zones(l)%speciation_alg%num_var_act_species) !> write component matrix elements [-]
                end do
                !write(unit,"(/,10x,'Initial concentration of components:'/)")  !> write initial component concentrations header
                !do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species !> loop over primary species
                !    write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:), &
                !        this%waters_init(wat_indices(j))%get_conc_var_act()), j=1,size(wat_indices)) !> write initial component concentration [C] computed from component matrix and species concentrations
                !end do
                !write(unit,"(/,10x,'Final concentration of components:'/)")  !> write final component concentrations header
                !do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species !> loop over primary species
                !    write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:), &
                !        this%waters(wat_indices(j))%get_conc_var_act()), j=1,size(wat_indices)) !> write final component concentration [C] computed from component matrix and species concentrations
                !end do
            end if
            ! write(unit,"(/,10x,'Final concentration of aqueous species:'/)")
            ! do i=1,this%chem_syst%aq_phase%num_species
            !     write(unit,"(10x,*(ES15.5))") (this%waters(wat_indices(j))%concentrations(this%waters(&
            !         wat_indices(j))%indices_aq_species(i)), j=1,size(wat_indices))
            ! end do
            ! write(unit,"(/,10x,'Final activities aqueous species:'/)")
            ! do i=1,this%chem_syst%aq_phase%num_species
            !     write(unit,"(10x,*(ES15.5))") (this%waters(wat_indices(j))%activities(this%waters(&
            !         wat_indices(j))%indices_aq_species(i)), j=1,size(wat_indices))
            ! end do
            !print *, this%target_solids_init(16)%mineral_zone%num_minerals_kin
            !print *, this%target_solids_init(16)%vol_fracts
            ! write(unit,"(/,10x,'Initial time of target solids:'/)")
            ! write(unit,"(10x,*(ES15.5))") (this%target_solids_init(tar_sol_indices(j))%time, j=1,size(tar_sol_indices))
            ! write(unit,"(/,10x,'Time of target solids:'/)")
            ! write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%time, j=1,size(tar_sol_indices))
            if (.not. mesh_too_large) then !> skip per-cell mineral/surface/gas concentrations & rates for large meshes (see binary mesh file)
            write(unit,"(/,10x,'Initial volumetric fractions of minerals:'/)")  !> write initial mineral volumetric fractions header
             do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in reactive zone
                 !print *,  this%target_solids_init(wat_indices(1))%vol_fracts
                 write(unit,"(10x,*(ES15.5))") (this%target_solids_init(tar_sol_indices(j))%min_pore_vol_ratios(this%target_solids_init(&
                      tar_sol_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(tar_sol_indices)) !> write initial mineral volumetric fraction φ_min [-]
             end do
             write(unit,"(/,10x,'Volumetric fractions of minerals:'/)")  !> write final mineral volumetric fractions header
             do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in reactive zone
                 write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%min_pore_vol_ratios(this%target_solids(&
                     tar_sol_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(tar_sol_indices)) !> write final mineral volumetric fraction φ_min [-]
             end do
             write(unit,"(/,10x,'Initial concentration of minerals:'/)")  !> write initial mineral concentrations header
             do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in reactive zone
                 write(unit,"(10x,*(ES15.5))") (this%target_solids_init(tar_sol_indices(j))%concentrations(this%target_solids_init(&
                     tar_sol_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(tar_sol_indices)) !> write initial mineral concentration [C]
             end do
             write(unit,"(/,10x,'Concentration of minerals:'/)")  !> write final mineral concentrations header
             do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in reactive zone
                 write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%concentrations(this%target_solids(&
                     tar_sol_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(tar_sol_indices)) !> write final mineral concentration [C]
             end do
            write(unit,"(/,10x,'Initial concentration of surface complexes:'/)")  !> write initial surface complex concentrations header
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
                write(unit,"(10x,*(ES25.15))") (this%target_solids_init(tar_sol_indices(j))%concentrations(&
                    this%target_solids_init(tar_sol_indices(j))%mineral_zone%num_minerals+i), j=1,size(tar_sol_indices)) !> write initial surface complex concentration [C]
            end do
            write(unit,"(/,10x,'Concentration of surface complexes:'/)")  !> write final surface complex concentrations header
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
                write(unit,"(10x,*(ES25.15))") (this%target_solids(tar_sol_indices(j))%concentrations(&
                    this%target_solids(tar_sol_indices(j))%mineral_zone%num_minerals+i), j=1,size(tar_sol_indices)) !> write final surface complex concentration [C]
            end do
            write(unit,"(/,10x,'Initial activity of surface complexes:'/)")  !> write initial surface complex activities header
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
                write(unit,"(10x,*(ES25.15))") (this%target_solids_init(tar_sol_indices(j))%activities(&
                this%target_solids_init(tar_sol_indices(j))%mineral_zone%num_minerals+i), j=1,size(tar_sol_indices)) !> write initial surface complex activity [-]
            end do
            write(unit,"(/,10x,'Activity of surface complexes:'/)")  !> write final surface complex activities header
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl !> loop over surface complexes in zone
                write(unit,"(10x,*(ES25.15))") (this%target_solids(tar_sol_indices(j))%activities(&
                    this%target_solids(tar_sol_indices(j))%mineral_zone%num_minerals+i), j=1,size(tar_sol_indices)) !> write final surface complex activity [-]
            end do
            write(unit,"(/,10x,'Initial concentration of gases:'/)")  !> write initial gas concentrations header
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq !> loop over equilibrium gases in zone
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(tar_sol_indices(j))%concentrations(i), j=1,size(tar_sol_indices)) !> write initial gas concentration [C]
            end do
            write(unit,"(/,10x,'Concentration of gases:'/)")  !> write final gas concentrations header
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq !> loop over equilibrium gases in zone
                write(unit,"(10x,*(ES15.5))") (this%target_gases(tar_sol_indices(j))%concentrations(i), &
                j=1,size(tar_sol_indices)) !> write final gas concentration [C]
            end do
            write(unit,"(/,10x,'Initial partial pressures of gases:'/)")  !> write initial gas partial pressures header
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq !> loop over equilibrium gases in zone
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(tar_sol_indices(j))%activities(i), j=1,size(tar_sol_indices)) !> write initial partial pressure P_i [atm]
            end do
            write(unit,"(/,10x,'Partial pressures of gases:'/)")  !> write final gas partial pressures header
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq !> loop over equilibrium gases in zone
                write(unit,"(10x,*(ES15.5))") (this%target_gases(tar_sol_indices(j))%activities(i),j=1,size(tar_sol_indices)) !> write final partial pressure P_i [atm]
            end do
            write(unit,"(/,10x,'Initial volume of gas:'/)")  !> write initial gas volume header
            if (this%reactive_zones(l)%gas_phase%num_species>0) then !> check if gas phase exists
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(tar_sol_indices(j))%volume, j=1,size(tar_sol_indices)) !> write initial gas volume [L³]
            end if
            write(unit,"(/,10x,'Volume of gas:'/)")  !> write final gas volume header
            if (this%reactive_zones(l)%gas_phase%num_species>0) then !> check if gas phase exists
                write(unit,"(10x,*(ES15.5))") (this%target_gases(tar_sol_indices(j))%volume, j=1,size(tar_sol_indices)) !> write final gas volume [L³]
            end if
            
            write(unit,"(/,10x,'Mineral equilibrium reaction amounts:'/)")  !> write mineral equilibrium reaction amounts header
            do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in zone
                write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_amnt_eq(i), j=1,size(tar_sol_indices)) !> write mineral equilibrium reaction amount R_e [mol]
            end do
            write(unit,"(/,10x,'Mineral equilibrium reaction rates:'/)")  !> write mineral equilibrium reaction rates header
            do i=1,this%reactive_zones(l)%num_minerals !> loop over minerals in zone
               write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_rate_eq(i), &
                    j=1,size(tar_sol_indices))
            end do
            write(unit,"(/,10x,'Exchange equilibrium reaction rates:'/)")  !> write cation exchange reaction rates header
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_exch_cats !> loop over exchangeable cations
                write(unit,"(10x,*(ES15.5))") (this%target_solids(tar_sol_indices(j))%rxn_rate_eq(&
                this%reactive_zones(l)%num_minerals+i), j=1,size(tar_sol_indices)) !> write mean exchange reaction rate r_e_mean [1/T]
            end do
            write(unit,"(/,10x,'Gas equilibrium reaction rates:'/)")  !> write gas equilibrium reaction rates header
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq !> loop over equilibrium gases
                write(unit,"(10x,*(ES15.5))") (this%target_gases(tar_sol_indices(j))%rxn_rate_eq(i), j=1,size(tar_sol_indices)) !> write mean gas equilibrium rate r_e_mean [1/T]
            end do
            end if !> mesh_too_large: end of per-cell mineral/surface/gas concentrations & rates
            !deallocate(wat_indices,tar_sol_indices) !> deallocate index arrays
        end if
    end do
    close(unit)
end subroutine write_chemistry