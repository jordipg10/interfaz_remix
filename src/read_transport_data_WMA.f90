!> @file read_transport_data_WMA.f90
!> @brief Reads Water Mixing Approach (WMA) transport data from file
!>
!> @details
!> This subroutine reads the mixing ratios and mixing water indices from the file '_WMA_lambdas.dat'.
!> The Water Mixing Approach represents advective transport through mixing ratios between
!> neighboring water parcels in a 1D domain.
!>
!> The file structure contains four data sections:
!> - MIXING RATIOS CONCENTRATIONS: Lambdas for each target water (transport step, all waters)
!> - MIXING RATIOS REACTIONS: Lambdas for each target water (reaction step, target waters only)
!> - MIXING WATERS CONCENTRATIONS: Indices of waters contributing to concentration mixing (target + external waters)
!> - MIXING WATERS REACTIONS: Indices of waters contributing to reaction mixing (target waters only)
!>
!> Key assumptions:
!> - Upstream/downstream classification is deduced from water indices (index < target = upstream)
!> - Concentration mixing includes boundary waters; reaction mixing includes only target waters
!> - If MIXING RATIOS REACTIONS / MIXING WATERS REACTIONS sections are absent,
!>   reaction data defaults to concentration data
!>
!> File format ('_WMA_lambdas.dat'):
!> ```
!> MIXING RATIOS CONCENTRATIONS
!> <num_ratios_1> <lambda_1> <lambda_2> ... <lambda_n>
!> ...
!> 0  ! terminator
!> MIXING RATIOS REACTIONS
!> <num_ratios_1> <lambda_1> <lambda_2> ... <lambda_m>
!> ...
!> 0  ! terminator
!> MIXING WATERS CONCENTRATIONS
!> <target_index> <mix_water_1> <mix_water_2> ... <mix_water_n>
!> ...
!> MIXING WATERS REACTIONS
!> <target_index> <mix_water_1> <mix_water_2> ... <mix_water_m>
!> ...
!> end
!> ```
!>
!> @param[in,out] this Transport object (transport_1D_transient_c) - modified with mixing data
!> @param[in] dir Directory for data file
!> @param[in] root Root name of the data file (file will be root//'_WMA_lambdas.dat')
!>
!> @note Upstream/downstream is determined by comparing each water index to the target water index
!>
!> @see allocate_mix_conc_indices
!>
!> Algorithm overview:
!> 1. Open '_WMA_lambdas.dat' file
!> 2. First pass: Count target waters and read dimensions for concentration and reaction mixing ratios
!> 3. Allocate spatial discretization, mixing ratio arrays and mixing index arrays
!> 4. Second pass: Read mixing ratio values and mixing water indices for both sections
!> 5. Reorder mixing data: separate current, upstream, and downstream waters (for both sections)
!>
!> @author Jordi Petchamé-Guerrero
!> @date 2025
!>
subroutine read_transport_data_WMA(this,dir,root)
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c     !< 1D Eulerian homogeneous mesh class
    use transport_transient_m, only: transport_1D_transient_c !< Transient 1D transport class
    use arrays_m, only: real_array_c, int_array_c            !< Real and integer array classes

    implicit none                                            !< Require explicit variable declarations
    class(transport_1D_transient_c) :: this                   !< Transport object (modified with mixing data)
    character(len=*), intent(in) :: dir                       !< Directory path for the data file
    character(len=*), intent(in) :: root                      !< Root name of data file (prefix for '_WMA_lambdas.dat')

    character(len=256) :: label                               !< Section label read from file header lines
    integer(kind=4) :: i                                     !< Loop counter over target waters
    integer(kind=4) :: j                                     !< Loop counter over mixing ratios within a target
    integer(kind=4) :: tar_dim                               !< Number of mixing ratios for current target (read from file)
    integer(kind=4) :: unit                                  !< Fortran file unit number for I/O operations
    integer(kind=4) :: num_targets                           !< Total number of target waters counted in pass 1
    integer(kind=4), allocatable :: dims_conc(:)             !< Array of concentration mixing ratio dimensions per target
    integer(kind=4), allocatable :: dims_R(:)                !< Array of reaction mixing ratio dimensions per target
    type(real_array_c) :: mixing_ratios_conc                  !< Temporary concentration mixing ratios (before reordering)
    type(real_array_c) :: mixing_ratios_R                     !< Temporary reaction mixing ratios (before reordering)
    logical :: has_react_ratios                                !< Flag: .true. if MIXING RATIOS REACTIONS section was found
    logical :: has_react_waters                                !< Flag: .true. if MIXING WATERS REACTIONS section was found

    unit = 59                                                !< Set file unit number for I/O operations
    has_react_ratios = .false.                                !< Initialise: assume no reaction ratios section
    has_react_waters = .false.                                !< Initialise: assume no reaction waters section

    !> ================================================================
    !> Pass 1: Count target waters and read mixing ratio dimensions
    !> ================================================================
    open(unit, file=dir//root//'_WMA_lambdas.dat', status='old', action='read') !< Open WMA lambdas data file for reading
    do                                                       !< Loop through file sections until 'end' marker
        read(unit,*) label                                   !< Read section label from current line
        if (label .eq. 'end') then                           !< Check for end-of-file marker
            exit                                             !< Exit pass 1 loop
        else if (label .eq. 'MIXING RATIOS CONCENTRATIONS') then !< Check for concentration mixing ratios section
            num_targets = 0                                  !< Initialise target water counter to zero
            do                                               !< Loop to count targets and record their dimensions
                read(unit,*) tar_dim                          !< Read number of mixing ratios for this target
                if (tar_dim .eq. 0) exit                     !< Zero terminates the section
                num_targets = num_targets + 1                !< Increment target water counter
                if (allocated(dims_conc)) then               !< Check if dims_conc already has entries
                    dims_conc = [dims_conc, tar_dim]          !< Append new dimension via auto-reallocation
                else                                         !< First target: initialise the array
                    dims_conc = [tar_dim]                     !< Create dims_conc with single element
                end if                                       !< End allocated check
            end do                                           !< End target counting loop
        else if (label .eq. 'MIXING RATIOS REACTIONS') then  !< Check for reaction mixing ratios section
            has_react_ratios = .true.                         !< Mark that reaction ratios section is present
            allocate(dims_R(num_targets))                     !< Allocate reaction dimensions array for all targets
            do i = 1, num_targets                            !< Loop over each target water
                read(unit,*) dims_R(i)                        !< Read number of reaction mixing ratios for target i
            end do                                           !< End reaction dimensions reading loop
            read(unit,*) tar_dim                              !< Read section terminator (expected 0)
        end if                                               !< End section label check
    end do                                                   !< End pass 1 main loop

    !> Fallback: if no reaction ratios section, default to concentration dimensions
    if (.not. has_react_ratios) then                          !< Check if reaction ratios section was absent
        dims_R = dims_conc(1:num_targets)                    !< Copy concentration dimensions as reaction dimensions
    end if                                                   !< End fallback check

    !> Allocate spatial discretization
    allocate(mesh_1D_Euler_homog_c :: this%spatial_discr)     !< Polymorphic allocation of 1D Euler mesh
    call this%spatial_discr%set_Num_targets(num_targets)      !< Set number of target waters in the mesh

    !> Allocate mixing ratio arrays (concentration and reaction)
    call mixing_ratios_conc%allocate_array(num_targets)       !< Allocate temp concentration mixing ratio columns
    call mixing_ratios_R%allocate_array(num_targets)          !< Allocate temp reaction mixing ratio columns
    do i = 1, num_targets                                    !< Loop over each target water
        call mixing_ratios_conc%cols(i)%allocate_vector(dims_conc(i)) !< Allocate concentration vector with dimension from file
        call mixing_ratios_R%cols(i)%allocate_vector(dims_R(i))       !< Allocate reaction vector with dimension from file
    end do                                                   !< End allocation loop
    deallocate(dims_conc, dims_R)                             !< Free temporary dimension arrays (no longer needed)

    !> Allocate concentration index structure
    this%mixing_ratios_conc = mixing_ratios_conc              !< Copy temp structure to object (needed for allocate_mix_conc_indices)
    call this%allocate_mix_conc_indices()                     !< Allocate mix_conc_indices with dim+2 per target (for up/down counts)

    !> Allocate reaction index structure
    this%mixing_ratios_R = mixing_ratios_R                   !< Copy temp structure to object for index allocation
    call this%mix_react_indices%allocate_array(num_targets)   !< Allocate reaction index columns
    do i = 1, num_targets                                    !< Loop over each target water
        call this%mix_react_indices%cols(i)%allocate_vector(mixing_ratios_R%cols(i)%dim + 2) !< Allocate with dim+2 for up/down counts
    end do                                                   !< End reaction index allocation loop

    rewind(unit)                                             !< Reset file position to beginning for pass 2

    !> ================================================================
    !> Pass 2: Read mixing ratio values and mixing water indices
    !> ================================================================
    do                                                       !< Loop through file sections until 'end' marker
        read(unit,*) label                                   !< Read section label from current line
        if (label .eq. 'end') then                           !< Check for end-of-file marker
            exit                                             !< Exit pass 2 loop
        else if (label .eq. 'MIXING RATIOS CONCENTRATIONS') then !< Check for concentration mixing ratios section
            do i = 1, num_targets                            !< Loop over each target water
                read(unit,*) tar_dim, (mixing_ratios_conc%cols(i)%col_1(j), j=1,tar_dim) !< Read dimension and lambda values
            end do                                           !< End concentration ratios reading loop
        else if (label .eq. 'MIXING RATIOS REACTIONS') then  !< Check for reaction mixing ratios section
            do i = 1, num_targets                            !< Loop over each target water
                read(unit,*) tar_dim, (mixing_ratios_R%cols(i)%col_1(j), j=1,tar_dim) !< Read dimension and lambda values
            end do                                           !< End reaction ratios reading loop
        else if (label .eq. 'MIXING WATERS CONCENTRATIONS') then !< Check for concentration mixing waters section
            call read_and_reorder_waters(num_targets, mixing_ratios_conc, & !< Read indices and reorder concentration data
                this%mixing_ratios_conc, this%mix_conc_indices)              !< into [self, upstream..., downstream...] layout
        else if (label .eq. 'MIXING WATERS REACTIONS') then  !< Check for reaction mixing waters section
            has_react_waters = .true.                         !< Mark that reaction waters section is present
            call read_and_reorder_waters(num_targets, mixing_ratios_R, & !< Read indices and reorder reaction data
                this%mixing_ratios_R, this%mix_react_indices)              !< into [self, upstream..., downstream...] layout
        end if                                               !< End section label check
    end do                                                   !< End pass 2 main loop

    !> Fallback: if reaction sections absent, default to concentration data
    if (.not. has_react_ratios) then                          !< No reaction ratios section => copy all conc data
        do i = 1, num_targets                                !< Loop over each target water
            this%mixing_ratios_R%cols(i)%col_1 = this%mixing_ratios_conc%cols(i)%col_1 !< Copy reordered conc ratios to reaction ratios
            this%mix_react_indices%cols(i)%col_1 = this%mix_conc_indices%cols(i)%col_1 !< Copy reordered conc indices to reaction indices
        end do                                               !< End full fallback copy loop
    else if (.not. has_react_waters) then                     !< Reaction ratios present but waters section absent
        do i = 1, num_targets                                !< Loop over each target water
            this%mix_react_indices%cols(i)%col_1 = this%mix_conc_indices%cols(i)%col_1 !< Copy conc indices to reaction indices
        end do                                               !< End indices-only fallback copy loop
    end if                                                   !< End fallback check

    close(unit)                                              !< Close the data file

contains

    !> @brief Reads mixing water indices and reorders ratios into [self, upstream..., downstream...] layout
    !>
    !> @details For each target water, reads the target index and its mixing water indices from
    !>          the file unit. Classifies each mixing water as upstream (index < target) or
    !>          downstream (index > target), then reorders both the mixing ratios and indices
    !>          into the canonical layout: [self, upstream_1..upstream_n, downstream_1..downstream_m].
    !>          The last two positions store the upstream and downstream counts respectively.
    !>
    !> @param[in]     n                Number of target waters to read
    !> @param[in]     raw_ratios       Mixing ratios as read from file (before reordering)
    !> @param[in,out] reordered_ratios Mixing ratios reordered into canonical layout
    !> @param[in,out] indices          Mixing water index structure reordered into canonical layout
    !>
    !> @note Accesses the parent subroutine's file unit via host association
    subroutine read_and_reorder_waters(n, raw_ratios, reordered_ratios, indices)
        integer(kind=4), intent(in) :: n                     !< Number of target waters to process
        type(real_array_c), intent(in) :: raw_ratios         !< Mixing ratios before reordering (from pass 2 read)
        type(real_array_c), intent(inout) :: reordered_ratios !< Mixing ratios after reordering (stored in object)
        type(int_array_c), intent(inout) :: indices          !< Mixing water index structure (stored in object)

        integer(kind=4) :: ii                                !< Loop counter over target waters
        integer(kind=4) :: jj                                !< Loop counter over mixing waters within a target
        integer(kind=4) :: mix_wat_ind                       !< Target water index read from file
        integer(kind=4) :: dim_i                             !< Dimension of mixing ratios for current target
        integer(kind=4) :: num_up                            !< Count of upstream mixing waters (index < target)
        integer(kind=4) :: num_down                          !< Count of downstream mixing waters (index > target)
        integer(kind=4) :: up_pos                            !< Placement counter for upstream entries during reorder
        integer(kind=4) :: down_pos                          !< Placement counter for downstream entries during reorder
        integer(kind=4), allocatable :: mix_wat_idx(:)       !< Temporary array of mixing water indices for current target

        do ii = 1, n                                         !< Loop over each target water
            dim_i = raw_ratios%cols(ii)%dim                  !< Get number of mixing waters for this target
            allocate(mix_wat_idx(dim_i))                     !< Allocate temp array for mixing water indices
            read(unit,*) mix_wat_ind, mix_wat_idx            !< Read target index and all mixing water indices from file
            indices%cols(ii)%col_1(1) = mix_wat_ind          !< Store self (target) index at position 1

            !> Count upstream and downstream waters by comparing index to target
            num_up = 0                                       !< Initialise upstream counter to zero
            num_down = 0                                     !< Initialise downstream counter to zero
            do jj = 1, dim_i                                 !< Loop over each mixing water
                if (mix_wat_idx(jj) < mix_wat_ind) then      !< Water index < target index => upstream
                    num_up = num_up + 1                      !< Increment upstream count
                else if (mix_wat_idx(jj) > mix_wat_ind) then !< Water index > target index => downstream
                    num_down = num_down + 1                  !< Increment downstream count
                end if                                       !< Equal to target = self, skip counting
            end do                                           !< End upstream/downstream counting loop
            indices%cols(ii)%col_1(indices%cols(ii)%dim - 1) = num_up  !< Store upstream count at second-to-last position
            indices%cols(ii)%col_1(indices%cols(ii)%dim) = num_down     !< Store downstream count at last position

            !> Reorder: self at position 1, upstream at 2..1+num_up, downstream at 2+num_up..
            up_pos = 0                                       !< Initialise upstream placement counter
            down_pos = 0                                     !< Initialise downstream placement counter
            do jj = 1, dim_i                                 !< Loop over each mixing water for reordering
                if (mix_wat_idx(jj) .eq. mix_wat_ind) then   !< Current water is self (target)
                    reordered_ratios%cols(ii)%col_1(1) = raw_ratios%cols(ii)%col_1(jj) !< Place self ratio at position 1
                else if (mix_wat_idx(jj) < mix_wat_ind) then !< Current water is upstream
                    up_pos = up_pos + 1                      !< Advance upstream placement position
                    reordered_ratios%cols(ii)%col_1(1 + up_pos) = raw_ratios%cols(ii)%col_1(jj) !< Place upstream ratio
                    indices%cols(ii)%col_1(1 + up_pos) = mix_wat_idx(jj)                        !< Place upstream index
                else                                         !< Current water is downstream
                    down_pos = down_pos + 1                  !< Advance downstream placement position
                    reordered_ratios%cols(ii)%col_1(1 + num_up + down_pos) = raw_ratios%cols(ii)%col_1(jj) !< Place downstream ratio
                    indices%cols(ii)%col_1(1 + num_up + down_pos) = mix_wat_idx(jj)                        !< Place downstream index
                end if                                       !< End upstream/downstream/self check
            end do                                           !< End reordering loop
            deallocate(mix_wat_idx)                           !< Free temporary mixing water indices array
        end do                                               !< End target waters loop
    end subroutine                                           !< End read_and_reorder_waters

end subroutine                                               !< End read_transport_data_WMA
