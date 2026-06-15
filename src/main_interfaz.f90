!> @file main_interfaz.f90
!> @brief Interactive driver for 1D reactive mixing iterations using the WMA
!>        (Water Mixing Approach), with explicit time integration of kinetics
!>        when applicable.
!> @details Reads chemistry/database paths and problem inputs from the user,
!>          loads the chemical system, writes initial/external water type
!>          component concentrations, and runs a loop of reactive mixing
!>          iterations using either a constant or variable time step.
program main_interfaz
    !> Import chemistry container type and the five reactive-mixing interface procedures:
    !>   - interfaz_esp_arch              : kinetic-only (no equilibrium reactions).
    !>   - interfaz_comps_arch_eq         : equilibrium-only, u_tilde with rows=components, cols=targets.
    !>   - interfaz_comps_arch_eq_T       : equilibrium-only, u_tilde transposed (rows=targets, cols=components).
    !>   - interfaz_comps_arch_eq_kin     : equilibrium + kinetic, u_tilde with rows=components, cols=targets.
    !>   - interfaz_comps_arch_eq_kin_T   : equilibrium + kinetic, u_tilde transposed (rows=targets, cols=components).
    use chemistry_m, only: chemistry_c, interfaz_comps_arch_eq_kin, interfaz_comps_arch_eq, interfaz_esp_arch, &
        interfaz_comps_arch_eq_kin_T, interfaz_comps_arch_eq_T, interfaz_comps_arch_eq_kin_mix
    !> Standard intrinsic module providing IEEE arithmetic helpers.
    use, intrinsic :: ieee_arithmetic
    !> Standard intrinsic module providing IEEE floating-point exception flags.
    use, intrinsic :: ieee_exceptions
    !> Disable implicit typing: every variable must be declared explicitly.
    implicit none
    !> Objects
    type(chemistry_c) :: my_chem !>< Chemistry container holding chemical system, waters and reactive zones.
    !> Variables
    integer(kind=4) :: num_aq_comps !>< Number of aqueous components in the chemical system.
    integer(kind=4) :: num_tar !>< Number of targets in the mesh.
    integer(kind=4) :: flag_transpose !>< Whether the input file has rows=targets & columns=components (1) or rows=components & columns=targets (0).
    integer(kind=4) :: flag_mix !>< Whether the input concentrations are PRIOR to mixing so the program performs the conservative mixing (1) or already after conservative transport (0).
    integer(kind=4) :: flag_wat_types !>< Whether file_u_wat_types has already been generated (1) or must be generated now (0).
    real(kind=8) :: Delta_t !>< Time step value (single reactive-mixing iteration).
    character(len=100) :: path_DB,path_pb,root_files,file_u_tilde,file_u_new,file_u_wat_types,file_mix !>< Raw fixed-length path and file name inputs.
    character(len=:), allocatable :: path_DB_trimmed, path_pb_trimmed, root_files_trimmed, file_u_tilde_trimmed, file_u_new_trimmed,&
        file_u_wat_types_trimmed, file_mix_trimmed !>< Trimmed deferred-length versions of the input strings.
    integer :: ios !>< I/O status code returned by safe read statements.
    logical :: has_ieee !>< True if the runtime supports the IEEE intrinsic standard.
    character(len=512) :: buf !>< Line buffer used by read_clean_line to skip comments/blank lines.
    procedure(interfaz_esp_arch), pointer :: p_interfaz=>null() !>< Procedure pointer dispatching to the appropriate reactive-mixing interface.
    !> Pre-Process
    !> Query whether IEEE intrinsic support is available in this build.
    has_ieee = ieee_support_standard()
    !> Configure IEEE behaviour and clear all relevant exception flags.
    if (has_ieee) then
        !> Disable gradual underflow (use abrupt underflow to zero).
        call ieee_set_underflow_mode(gradual=.false.)
        !> Clear underflow exception flag.
        call ieee_set_flag(ieee_underflow, .false.)
        !> Clear inexact-result exception flag.
        call ieee_set_flag(ieee_inexact, .false.)
        !> Clear overflow exception flag.
        call ieee_set_flag(ieee_overflow, .false.)
        !> Clear divide-by-zero exception flag.
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        !> Clear invalid-operation exception flag.
        call ieee_set_flag(ieee_invalid, .false.)
    end if
    !> Process
    !> Banner describing what the interactive program does.
    write(*,*) "This is the interface to solve the reactive part of a 1D reactive mixing iteration using &
        explicit Euler."
    !> Prompt the user for the database path.
    write(*,*) "Database path: "
    !> Read the database path from standard input (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) path_DB
    !> If the read failed (e.g. EOF when running non-interactively), abort gracefully.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading path_DB. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the database path string.
    path_DB_trimmed = trim(path_DB)
    !> Prompt the user for the problem path.
    write(*,*) "Problem path: "
    !> Read the problem path from standard input (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) path_pb
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading path_pb. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the problem path string.
    path_pb_trimmed = trim(path_pb)
    !> Prompt the user for the common root of input/output filenames.
    write(*,*) "Root of the input and output files: "
    !> Read the file root from standard input (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) root_files
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading root_files. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the file root string.
    root_files_trimmed = trim(root_files)
    !> Prompt the user for the number of targets in the mesh.
    write(*,*) "How many targets does the mesh have?"
    !> Read the number of targets (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) num_tar
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading num_tar. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    end if
    !> we read chemistry
    !> Load the chemical system, water types and reactive zones from disk.
    !> num_tar is passed so that the target waters file can be validated against the mesh size.
    call my_chem%read_chemistry_interface(root_files_trimmed,path_pb_trimmed,path_DB_trimmed,num_tar)
    !> Ask whether the water-types concentration file has already been generated
    !> in a previous run. If yes (1) we just read its name to display it; if no
    !> (0) we generate it now via write_conc_comp_tar_wat.
    write(*,*) "Has the file with the aqueous component concentrations of the initial and external water types &
        already been generated? (1: yes, 0: no)"
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) flag_wat_types
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading flag_wat_types. Run in an interactive terminal or redirect from fort.5'; call &
            safe_stop(1)
    else if (flag_wat_types /= 0 .and. flag_wat_types /= 1) then
        error stop "Invalid option. It must be 1 (file already generated) or 0 (generate it now)."
    end if
    !> Prompt for the file name (used either as output to write, or just as info).
    if (flag_wat_types == 0) then
        write(*,*) "Name of the file where you want me to write the aqueous component concentrations &
            of the initial and external water types?"
    else
        write(*,*) "Name of the already-generated file with the water-type concentrations:"
    end if
    !> Read the water-types filename (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) file_u_wat_types !> file with initial and external water types
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading file_u_wat_types. Run in an interactive terminal or redirect from fort.5'; call & 
            safe_stop(1)
    end if
    !> Trim trailing blanks from the water-types output filename.
    file_u_wat_types_trimmed = trim(file_u_wat_types)
    !> Query the first target water for its number of aqueous components.
    !> This value is passed to the reactive-mixing interface so it can correctly
    !> read the u_tilde input file (one entry per aqueous component per target).
    num_aq_comps=my_chem%get_num_aq_comps_tar_wat()
    !> Echo the number of aqueous components for diagnostics.
    !print *, "Number of aqueous components: ", num_aq_comps
    !> Generate the water-types concentration file only when requested.
    if (flag_wat_types == 0) then
        !> we write concentrations of target waters (domain and external)
        call my_chem%write_conc_comp_tar_wat(path_pb_trimmed,file_u_wat_types_trimmed)
        write(*,*) 'File ' // trim(file_u_wat_types_trimmed) // ' generated successfully.'
    else
        write(*,*) 'Using the existing file ' // trim(file_u_wat_types_trimmed) // '.'
    end if
    !> write file name for u_tilde and u_new
    !> Prompt for the output file name for post-reactive-mixing concentrations.
    write(*,*) "Name of the file where you want me to write the concentrations after the reactive mixing?"
    !> Read the output filename for post-mixing concentrations (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) file_u_new
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading file_u_new. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the post-mixing output filename.
    file_u_new_trimmed = trim(file_u_new)
    !> Prompt for the input file containing u_tilde (post conservative-transport concentrations).
    write(*,*) "Name of the file containing the aqueous component concentrations after solving one iteration of &
        conservative transport? &
        IMPORTANT: The file must be located in the problem path."
    !> Read the input filename providing u_tilde (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) file_u_tilde !> file with u_tilde
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading file_u_tilde. Run in an interactive terminal or redirect from fort.5'; call &
            safe_stop(1)
    end if
    !> Ask whether the input concentrations are PRIOR to mixing. If so, the program
    !> performs the conservative mixing internally using a mixing-ratios file
    !> (rows=components, columns=waters is the fixed layout in this case, so the
    !> transpose question below is skipped).
    write(*,*) "Are the concentrations in " // trim(file_u_tilde) // " PRIOR to mixing, so I should also perform the &
        conservative mixing using a mixing ratios file? (1: yes, 0: no, they are already after conservative transport)"
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) flag_mix
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading flag_mix. Run in an interactive terminal or redirect from fort.5'; call &
            safe_stop(1)
    else if (flag_mix /= 0 .and. flag_mix /= 1) then
        error stop "Invalid option. It must be 1 (concentrations prior to mixing) or 0 (after conservative transport)."
    end if
    if (flag_mix == 1) then
        !> Prompt for the mixing-ratios file (row-major order: rows=target waters, columns=waters).
        write(*,*) "Name of the file containing the mixing ratios in row-major order (rows: target waters, &
            columns: waters)? IMPORTANT: The file must be located in the problem path and its column ordering &
            must match the columns (waters) of " // trim(file_u_tilde) // "."
        call read_clean_line(buf, ios)
        if (ios == 0) read(buf, *, iostat=ios) file_mix
        if (ios /= 0) then
            write(*,*) 'Error/EOF reading file_mix. Run in an interactive terminal or redirect from fort.5'; call &
                safe_stop(1)
        end if
        !> Layout is fixed for the mixing variant (rows=components, columns=waters).
        flag_transpose = 0
    else
        !> Ask the user whether the input matrix is transposed (rows=targets, columns=components).
        write(*,*) "Does the file " // trim(file_u_tilde) // " have rows as targets and columns as components? (1: yes, 0: no)"
        call read_clean_line(buf, ios)
        if (ios == 0) read(buf, *, iostat=ios) flag_transpose
        if (ios /= 0) then
            write(*,*) 'Error/EOF reading flag_transpose. Run in an interactive terminal or redirect from fort.5'; call &
                safe_stop(1)
        else if (flag_transpose /= 0 .and. flag_transpose /= 1) then
            error stop "Invalid option. It must be 1 (transposed component matrix) or 0 (not transposed)."
        end if
    end if
    !> Prompt for the time step value (single iteration).
    write(*,*) "Time step?"
    !> Read the time step (skipping comment/blank lines).
    call read_clean_line(buf, ios)
    if (ios == 0) read(buf, *, iostat=ios) Delta_t !> time step
    !> Validate the read: report I/O failure or non-positive time step.
    if (ios /= 0) then
        write(*,*) 'Error/EOF reading Delta_t. Run in an interactive terminal or redirect from fort.5'; call safe_stop(1)
    else if (Delta_t <= 0d0) then
        !> Time step must be strictly positive.
        error stop "The time step must be greater than zero"
    end if
    !> Trim trailing blanks from the u_tilde input filename.
    file_u_tilde_trimmed = trim(file_u_tilde) !> we trim file name
    !> When the concentrations are prior to mixing, the program performs the
    !> conservative mixing internally via interfaz_comps_arch_eq_kin_mix, which
    !> takes an extra mixing-ratios file argument and therefore cannot share the
    !> p_interfaz procedure pointer. Dispatch it directly and skip the pointer
    !> selection below.
    if (flag_mix == 1) then
        file_mix_trimmed = trim(file_mix)
        associate(tw0 => my_chem%waters(my_chem%tar_wat_indices(1)))
            if (tw0%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions == 0) then
                error stop "The conservative+reactive mixing variant requires equilibrium reactions."
            end if
        end associate
        write(*,*) "Proceeding with one conservative + reactive mixing iteration."
        call interfaz_comps_arch_eq_kin_mix(my_chem,path_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,&
            file_mix_trimmed,Delta_t,file_u_new_trimmed)
        !> Post-Process
        if (has_ieee) call clear_ieee_flags()
        write(*,*) "The program has finished."
        call safe_stop(0)
    end if
    !> Select the reactive-mixing interface via procedure pointer based on the
    !> chemical system and on the orientation of the u_tilde input file:
    !>   - no equilibrium reactions                  → interfaz_esp_arch (kinetic only)
    !>   - equilibrium reactions, no kinetics:
    !>       * flag_transpose == 0 (rows=components) → interfaz_comps_arch_eq
    !>       * flag_transpose == 1 (rows=targets)    → interfaz_comps_arch_eq_T
    !>   - equilibrium reactions AND kinetics:
    !>       * flag_transpose == 0 (rows=components) → interfaz_comps_arch_eq_kin
    !>       * flag_transpose == 1 (rows=targets)    → interfaz_comps_arch_eq_kin_T
    associate(tw0 => my_chem%waters(my_chem%tar_wat_indices(1)))
        if (tw0%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions == 0) then
            p_interfaz => interfaz_esp_arch
        else if (tw0%solid_chemistry%reactive_zone%chem_syst%num_kin_reacts == 0) then
            if (flag_transpose == 1) then
                p_interfaz => interfaz_comps_arch_eq_T
            else
                p_interfaz => interfaz_comps_arch_eq
            end if
        else
            if (flag_transpose == 1) then
                p_interfaz => interfaz_comps_arch_eq_kin_T
            else
                p_interfaz => interfaz_comps_arch_eq_kin
            end if
        end if
    end associate
    !> Run a single reactive-mixing iteration using the selected interface.
    write(*,*) "Proceeding with one reactive mixing iteration."
    call p_interfaz(my_chem,path_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
    !> Post-Process
    !> Clear any pending IEEE flags before normal termination.
    if (has_ieee) call clear_ieee_flags()
    !> Notify the user that execution finished successfully.
    write(*,*) "The program has finished."
contains
    !> @brief Read the next non-blank, non-comment line from standard input.
    !> @details Skips lines that are empty (after trimming) or whose first non-blank
    !>          character is '!' or '#'. This lets the user annotate input files
    !>          (e.g. fort.5 redirection) with comments without breaking the
    !>          list-directed reads in the main program.
    !>          Only FULL-LINE comments are stripped; '!' or '#' that appear in
    !>          the middle of a line are kept verbatim, so file paths containing
    !>          those characters are preserved. Place comments on their own line.
    !>          The caller is expected to parse @p line with an internal read,
    !>          e.g.  read(line, *, iostat=ios) value.
    !> @param[out] line Buffer that receives the next clean input line.
    !> @param[out] ios  I/O status: 0 on success, /= 0 on read error or end of file.
    subroutine read_clean_line(line, ios)
        character(len=*), intent(out) :: line
        integer, intent(out)          :: ios
        character(len=:), allocatable :: trimmed
        integer :: n
        do
            read(*, '(A)', iostat=ios) line
            if (ios /= 0) return                 !> propagate error / EOF to the caller
            trimmed = adjustl(line)
            n = len_trim(trimmed)
            if (n == 0) cycle                    !> skip blank lines
            if (trimmed(1:1) == '!' .or. trimmed(1:1) == '#') cycle  !> skip full-line comments
            line = trimmed
            return
        end do
    end subroutine read_clean_line

    !> @brief Clear all IEEE floating-point exception flags if IEEE support is available.
    !> @details Used both before normal termination and inside @ref safe_stop to avoid
    !>          spurious end-of-program warnings about pending floating-point exceptions.
    subroutine clear_ieee_flags()
        !> Only act if the IEEE intrinsic standard is supported by the runtime.
        if (ieee_support_standard()) then
            !> Clear underflow flag.
            call ieee_set_flag(ieee_underflow, .false.)
            !> Clear inexact flag.
            call ieee_set_flag(ieee_inexact, .false.)
            !> Clear overflow flag.
            call ieee_set_flag(ieee_overflow, .false.)
            !> Clear divide-by-zero flag.
            call ieee_set_flag(ieee_divide_by_zero, .false.)
            !> Clear invalid-operation flag.
            call ieee_set_flag(ieee_invalid, .false.)
        end if
    end subroutine clear_ieee_flags

    !> @brief Stop the program with a given exit code, first clearing IEEE flags.
    !> @param[in] code Integer exit/stop code passed to the Fortran @c stop statement.
    subroutine safe_stop(code)
        !> Exit code propagated to the operating system.
        integer, intent(in) :: code
        !> Clear pending IEEE flags so the runtime does not emit warnings.
        call clear_ieee_flags()
        !> Terminate execution with the requested code.
        stop code
    end subroutine safe_stop
end program main_interfaz