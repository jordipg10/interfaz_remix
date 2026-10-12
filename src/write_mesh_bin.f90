!> \file write_mesh_bin.f90
!> \brief Write per-cell mesh concentrations and reaction rates to a compact binary file.
!> \details
!>   For large meshes (e.g. the 8000-cell Paulina grid) the formatted *_chemistry.out file
!>   is unwieldy. This routine writes the per-target-water (per-cell) data needed for plotting
!>   to a raw little-endian float64 stream file, plus a small text metadata sidecar that
!>   records the dimensions, the species/mineral names and the block layout. The binary file
!>   contains ONLY float64 values, so it can be read back trivially, e.g. with numpy:
!>
!>     import numpy as np
!>     # parse <root>_..._mesh_meta.txt for: Nx Ny n_tw n_sp n_eq n_aqkin n_min
!>     d = np.fromfile('<root>_..._mesh.bin', dtype='<f8')
!>     o = 0
!>     x          = d[o:o+n_tw]; o += n_tw
!>     y          = d[o:o+n_tw]; o += n_tw
!>     conc       = d[o:o+n_sp*n_tw].reshape(n_sp, n_tw);       o += n_sp*n_tw
!>     aq_eq_rate = d[o:o+n_eq*n_tw].reshape(n_eq, n_tw);       o += n_eq*n_tw
!>     aq_kin_rate= d[o:o+n_aqkin*n_tw].reshape(n_aqkin, n_tw); o += n_aqkin*n_tw
!>     min_rate   = d[o:o+n_min*n_tw].reshape(n_min, n_tw);     o += n_min*n_tw
!>
!>   Every block is row = quantity index, column = target water (cell) in tar_wat_indices order
!>   (row-major mesh ordering: cell index = (row-1)*Nx + col). For a 2D mesh use x,y (or the
!>   Nx,Ny dims) to lay the cells out on the grid.
!>
!> \param[in] this RT_c reactive transport object (after the simulation has run)
!> \param[in] dir  directory for output files
!> \param[in] root root name for output files
!> \author Jordi Petchame-Guerrero
!> \date June 2026
subroutine write_mesh_bin(this, dir, root)
    use RT_m, only: RT_c, RT_1D_transient_c, RT_2D_transient_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c
    implicit none
    class(RT_c), intent(in) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root

    integer(kind=4) :: i, j, p, g, ndim
    integer(kind=4) :: n_tw, n_sp, n_eq, n_aqkin, n_min, Nx, Ny
    integer(kind=4) :: nmk, ntot
    integer(kind=4) :: ubin, umeta
    integer(kind=4) :: num_time                !< total number of time steps in the run (appended to filename)
    logical :: lagr_flag                       !< Lagrangian (vs Eulerian) transport flag (from the transient transport object)
    character(len=:), allocatable :: wma_tag, tpt_tag, nt_tag, base
    character(len=32) :: nt_buf                !< scratch buffer for formatting the time-step count
    real(kind=8), allocatable :: val(:)        !< per-species/rate work row over all cells
    real(kind=8), allocatable :: min_rate(:,:) !< (n_min, n_tw) mineral rates mapped to global minerals
    real(kind=8), allocatable :: xc(:), yc(:)  !< cell coordinates

    !> Dimensions (from chemistry/chemical system; available on the RT_c base type)
    n_tw    = this%chemistry%num_target_waters
    n_sp    = this%chemistry%chem_syst%aq_phase%num_species
    n_eq    = this%chemistry%chem_syst%aq_phase%num_aq_complexes
    n_aqkin = this%chemistry%chem_syst%num_aq_kin_reacts
    n_min   = this%chemistry%chem_syst%num_minerals

    !> Transport-dependent quantities (Lagr_flag, grid dims, cell coordinates) live on the
    !> concrete transient types, not on the abstract RT_c base, so reach them via select type.
    allocate(xc(n_tw), yc(n_tw))
    xc = 0d0
    yc = 0d0
    Nx = n_tw
    Ny = 1
    ndim = 1
    lagr_flag = .false.
    num_time = 0
    select type (this)
    type is (RT_1D_transient_c)
        lagr_flag = this%transport%Lagr_flag
        num_time = this%transport%time_discr%Num_time
        call gather_coords(this%transport%spatial_discr)
    type is (RT_2D_transient_c)
        lagr_flag = this%transport%Lagr_flag
        num_time = this%transport%time_discr%Num_time
        select type (mesh => this%transport%spatial_discr)
        type is (mesh_2D_Euler_homog_c)
            Nx = mesh%Num_cells_x
            Ny = mesh%Num_cells_y
            ndim = 2
        end select
        call gather_coords(this%transport%spatial_discr)
    class default
        error stop "write_mesh_bin: unsupported RT type (expected a transient RT object)"
    end select

    !> Output-name tags consistent with write_chemistry
    if (this%chemistry%lump_flag) then
        wma_tag = '_lumped'
    else
        wma_tag = '_consistent'
    end if
    if (lagr_flag) then
        tpt_tag = '_Lagr'
    else
        tpt_tag = '_Euler'
    end if
    if (num_time > 0) then
        write(nt_buf,'(I0)') num_time
        nt_tag = '_'//trim(nt_buf)//'steps'
    else
        nt_tag = ''
    end if
    base = dir//root//tpt_tag//wma_tag//nt_tag//'_mesh'

    !> --- Metadata sidecar (text) ---
    open(newunit=umeta, file=base//'_meta.txt', status='unknown', form='formatted')
    write(umeta,'(A)') 'REMIX_MESH_BIN 1'
    write(umeta,'(A)') '# dims: Nx Ny n_tw n_sp n_eq n_aqkin n_min'
    write(umeta,'(7(I0,1x))') Nx, Ny, n_tw, n_sp, n_eq, n_aqkin, n_min
    write(umeta,'(A)') 'SPECIES'
    do i = 1, n_sp
        write(umeta,'(A)') trim(this%chemistry%chem_syst%aq_phase%aq_species(i)%name)
    end do
    write(umeta,'(A)') 'MINERALS'
    do i = 1, n_min
        write(umeta,'(A)') trim(this%chemistry%chem_syst%minerals(i)%name)
    end do
    write(umeta,'(A)') 'LAYOUT  (name dtype rows cols ; binary order, float64 little-endian)'
    write(umeta,'(A,I0)')        'x           f8 1 ',       n_tw
    write(umeta,'(A,I0)')        'y           f8 1 ',       n_tw
    write(umeta,'(A,I0,1x,I0)')  'conc        f8 ', n_sp,    n_tw
    write(umeta,'(A,I0,1x,I0)')  'aq_eq_rate  f8 ', n_eq,    n_tw
    write(umeta,'(A,I0,1x,I0)')  'aq_kin_rate f8 ', n_aqkin, n_tw
    write(umeta,'(A,I0,1x,I0)')  'min_rate    f8 ', n_min,   n_tw
    write(umeta,'(A)') '# rows = quantity index, cols = cell (tar_wat order, cell=(row-1)*Nx+col)'
    write(umeta,'(A)') '# min_rate row m = global mineral m; 0 where that mineral is absent in a cell'
    close(umeta)

    !> --- Binary data (pure float64 stream) ---
    open(newunit=ubin, file=base//'.bin', status='unknown', form='unformatted', access='stream')

    !> Coordinates
    write(ubin) xc
    write(ubin) yc

    allocate(val(n_tw))

    !> Concentrations: one row per aqueous species, over all cells
    do i = 1, n_sp
        do j = 1, n_tw
            val(j) = this%target_waters(j)%aq_chem%concentrations( &
                this%target_waters(j)%aq_chem%indices_aq_species(i))
        end do
        write(ubin) val
    end do

    !> Aqueous equilibrium reaction rates (re_mean), one row per aqueous complex
    do i = 1, n_eq
        do j = 1, n_tw
            if (allocated(this%target_waters(j)%aq_chem%rxn_rate_eq) .and. &
                i <= size(this%target_waters(j)%aq_chem%rxn_rate_eq)) then
                val(j) = this%target_waters(j)%aq_chem%rxn_rate_eq(i)
            else
                val(j) = 0d0
            end if
        end do
        write(ubin) val
    end do

    !> Aqueous kinetic reaction rates (rk_mean), one row per aqueous kinetic reaction
    do i = 1, n_aqkin
        do j = 1, n_tw
            if (allocated(this%target_waters(j)%aq_chem%rxn_rate_kin_avg) .and. &
                i <= size(this%target_waters(j)%aq_chem%rxn_rate_kin_avg)) then
                val(j) = this%target_waters(j)%aq_chem%rxn_rate_kin_avg(i)
            else
                val(j) = 0d0
            end if
        end do
        write(ubin) val
    end do

    !> Mineral reaction rates mapped to the global mineral list.
    !> Per cell, each local mineral (kinetic first, then equilibrium, per ind_min_chem_syst)
    !> is placed at its global mineral slot; absent minerals stay 0.
    allocate(min_rate(n_min, n_tw))
    min_rate = 0d0
    if (n_min > 0) then
        do j = 1, n_tw
            associate(sc => this%target_waters(j)%aq_chem%solid_chemistry)
            if (associated(sc)) then
                if (associated(sc%mineral_zone)) then
                    nmk  = sc%mineral_zone%num_minerals_kin
                    ntot = sc%mineral_zone%num_minerals
                    do p = 1, ntot
                        if (.not. allocated(sc%mineral_zone%ind_min_chem_syst)) exit
                        if (p > size(sc%mineral_zone%ind_min_chem_syst)) exit
                        g = sc%mineral_zone%ind_min_chem_syst(p)
                        if (g < 1 .or. g > n_min) cycle
                        if (p <= nmk) then
                            if (allocated(sc%rxn_rate_kin_avg) .and. p <= size(sc%rxn_rate_kin_avg)) &
                                min_rate(g, j) = sc%rxn_rate_kin_avg(p)
                        else
                            if (allocated(sc%rxn_rate_eq) .and. (p - nmk) <= size(sc%rxn_rate_eq)) &
                                min_rate(g, j) = sc%rxn_rate_eq(p - nmk)
                        end if
                    end do
                end if
            end if
            end associate
        end do
        do i = 1, n_min
            write(ubin) min_rate(i, :)
        end do
    end if

    close(ubin)

    write(*,'(A)') ' Wrote mesh binary: '//base//'.bin (+ _meta.txt)'

    deallocate(val, xc, yc, min_rate)

contains

    !> Fill the host xc/yc arrays from the mesh target coordinates (row-major cell order).
    subroutine gather_coords(sd)
        use spatial_discr_m, only: spatial_discr_c
        class(spatial_discr_c), intent(in) :: sd
        integer(kind=4) :: jj
        do jj = 1, n_tw
            if (allocated(sd%targets(jj)%coord)) then
                xc(jj) = sd%targets(jj)%coord(1)
                if (ndim >= 2 .and. size(sd%targets(jj)%coord) >= 2) then
                    yc(jj) = sd%targets(jj)%coord(2)
                end if
            end if
        end do
    end subroutine gather_coords

end subroutine write_mesh_bin
