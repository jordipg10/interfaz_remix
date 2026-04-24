!> @file write_transport_data_WMA.f90
!> @brief Writes transport data from 1D reactive transport problem solved with Water Mixing Approach (WMA)
!> @details This subroutine outputs key transport simulation parameters to an already-opened file.
!> The Water Mixing Approach (WMA) is a Lagrangian-Eulerian method for reactive transport that
!> tracks mixing ratios of water parcels to solve advection-dominated problems efficiently.
!> This routine writes spatial discretization (target coordinates), temporal discretization
!> (time step, final time), and optionally mixing ratio and Lagrangian particle information.
!>
!> @par Output Format:
!> The output file contains:
!>   - Number of spatial targets (observation points)
!>   - Target coordinates [L]
!>   - Time step size Δt [T]
!>   - Final simulation time [T]
!>   - (Commented) Mixing ratios, reaction amounts, water indices, Lagrangian particle data
!>
!> @par Water Mixing Approach (WMA):
!> The WMA tracks mixing ratios of water parcels rather than concentrations directly, which
!> can be more efficient for advection-dominated reactive transport. Each cell maintains
!> mixing ratios indicating contributions from upstream water types.
!>
!> @see De Simoni et al. (2005), Water Resources Research, 41(11), W11410
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

!> @brief Write transport data for WMA-based 1D reactive transport simulation
!> @param[in] this 1D transient transport object containing spatial/temporal discretization and mixing ratios
!> @param[in] unit File unit number for output (file must be already opened by caller)
subroutine write_transport_data_WMA_2D(this,path,root)
    use transport_transient_m, only: transport_2D_transient_c !< Import 2D transient transport class definition
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c !< Import 2D homogeneous Eulerian mesh class definition
    implicit none !< Enforce explicit variable declarations for type safety
    class(transport_2D_transient_c), intent(in) :: this !< 2D transient transport object with WMA data (spatial_discr, time_discr, mixing_ratios) [-]
    character(len=*), intent(in) :: path !< Directory path for output files [-]
    character(len=*), intent(in) :: root !< Root name for output files [-]
    
    integer(kind=4) :: i !< Loop counter for iterating over targets or data rows [-]
    integer(kind=4) :: n !< Number of computational targets [-]
    integer(kind=4) :: unit !< File unit for results output [-]

    !> Compute number of computational cells accounting for boundary/target flags
    !> targets_flag indicates special targets (e.g., boundaries) that should be excluded from cell count
    n=this%spatial_discr%Num_targets !< Calculate active cells: total targets - special flags [-]
    
    !> Write transport data to transport file
    open(newunit=unit, file=path//root//'_transport.out', status='unknown', form='formatted')
    if (this%Lagr_flag) then
        write(unit,"(2x,'Method: Lagrangian'/)")
    else
        write(unit,"(2x,'Method: Eulerian'/)")
    end if
    write(unit,"(2x,'Type of targets:'/)")
    if (this%spatial_discr%targets_flag.eq.0) then
        write(unit,"(10x,'Cells')") !< Output target type
    else if (this%spatial_discr%targets_flag.eq.1) then
        write(unit,"(10x,'Edges')") !< Output target type
    else
        error stop 'Error in write_transport_data_WMA_2D: Unknown targets_flag value.'
    ! else if (this%spatial_discr%targets_flag.eq.2) then
    !     write(unit,"(2x,'Including source/sink targets')") !< Output target type
    ! else
    !     write(unit,"(2x,'Including boundary and source/sink targets')") !< Output target type
    end if
    !> Write number of spatial targets (observation points) to output file
    write(unit,"(2x,'Number of targets:',I5/)") n !< Output total number of targets (including boundaries)
    
    !> @note Legacy commented code: Alternative format writing individual target coordinates
    !do i=1,this%spatial_discr%Num_targets !< Would loop over all targets
        !write(unit,"(2x,'Target ',I5,': coord=',ES15.5)") i, this%spatial_discr%targets(i)%coord !< Would write: "Target i: coord=x" [L]
    !end do
    select type (mesh=>this%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        write(unit,"(2x,'Number of cells in x direction:',I5/)") mesh%Num_cells_x !< Output total number of targets (including boundaries)
        write(unit,"(2x,'Number of cells in y direction:',I5/)") mesh%Num_cells_y !< Output total number of targets (including boundaries)
        write(unit,"(2x,'Delta x:',ES15.5/)") mesh%Delta_x !< Output uniform cell size in x direction [L]
        write(unit,"(2x,'Delta y:',ES15.5/)") mesh%Delta_y !< Output uniform cell size in y direction [L]
    type is (mesh_2D_Euler_heterog_c)
        !write(unit,"(2x,'Number of cells in x direction:',I5/)") mesh%Num_cells_x !< Output total number of targets (including boundaries)
        !write(unit,"(2x,'Number of cells in y direction:',I5/)") mesh%Num_cells_y !< Output total number of targets (including boundaries)
    end select
    !> Write dispersion tensor for homogeneous flux case
    if (this%tpt_props_heterog%cst_flux_flag) then
        write(unit,"(2x,'Longitudinal dispersivity:',ES15.5/)") this%tpt_props_heterog%long_dispersivity
        write(unit,"(2x,'Transverse dispersivity:',ES15.5/)") this%tpt_props_heterog%transv_dispersivity
        if (allocated(this%tpt_props_heterog%disp_tensor_cent)) then
            write(unit,"(2x,'Dispersion tensor (homogeneous flux):'/)")
            write(unit,"(4x,'D_xx:',ES15.5)") this%tpt_props_heterog%disp_tensor_cent(1,1,1)
            write(unit,"(4x,'D_yy:',ES15.5)") this%tpt_props_heterog%disp_tensor_cent(1,1,2)
            write(unit,"(4x,'D_xy:',ES15.5/)") this%tpt_props_heterog%disp_tensor_cent(1,1,3)
        end if
    end if
    !> Write header for target coordinates section
    write(unit,"(2x,'Target coordinates (x,y):'/)")
    !> Write all target coordinates in a single line using implied-do loop
    !> Format: space, then all coordinates in scientific notation (ES15.5)
    do i=1,n !< Loop over all targets
        write(unit,"(2x,ES15.5,ES15.5)") this%spatial_discr%targets(i)%coord(1), &
            this%spatial_discr%targets(i)%coord(2) !< Output (x,y) coordinates [L] for all targets
    end do
    !> Write header for target coordinates section
    !write(unit,"(2x,'Target coordinates:'/)")  !< Output section header with blank line after
    
    !> Write all target coordinates in a single line using implied-do loop
    !> Format: space, then all coordinates in scientific notation (ES15.5)
    !write(unit,"(2x,*(ES15.5))") (this%spatial_discr%targets(i)%coord, i=1,this%spatial_discr%Num_targets) !< Output x-coordinates [L] for all targets
    
    !> Write header for time step section
    write(unit,"(/,2x,'Time step:'/)")  !< Output section header with blank lines before and after
    
    !> Write time step size retrieved from time discretization object
    write(unit,"(2x,ES15.5/)") this%time_discr%get_Delta_t() !< Output Δt [T] in scientific notation (may be constant or variable)
    
    !> Write header for final time section
    write(unit,"(2x,'Final time:'/)")  !< Output section header with blank line after
    
    !> Write final simulation time
    write(unit,"(2x,ES15.5/)") this%time_discr%Final_time !< Output t_final [T] in scientific notation
    
    !> @note Legacy commented code: Mixing ratios concentration data output
    !> @details Would output dimension and concentration mixing ratios for each column (cell/target)
    !> Format: dimension, then all mixing ratio values including boundary and source/sink terms
    ! write(unit,"(/,2x,'Dimension + Mixing ratios concentrations (by rows) (including boundary and sink/source terms):'/)")
    ! do i=1,this%mixing_ratios_conc%num_cols !< Loop over all columns (cells) in mixing ratios matrix
    !     write(unit,"(2x,I5,*(ES15.5))") this%mixing_ratios_conc%cols(i)%dim, (this%mixing_ratios_conc%cols(i)%col_1(j), &
    !         j=1,this%mixing_ratios_conc%cols(i)%dim) !< Output: dim, then mixing_ratio(1:dim) [-]
    ! end do
    
    !> @note Legacy commented code: Initial mixing ratios for reaction amounts
    !> @details Would output initial state of reaction amount mixing ratios at t=0
    ! write(unit,"(/,2x,'Initial mixing ratios reaction amounts (by rows):'/)")
    ! do i=1,this%mixing_ratios_R%num_cols !< Loop over all columns for initial reaction mixing ratios
    !     write(unit,"(2x,*(ES15.5))") (this%mixing_ratios_R_init%cols(i)%col_1(j), &
    !         j=1,this%mixing_ratios_R%cols(i)%dim) !< Output: initial R_k mixing ratios [M/L³ or dimensionless]
    ! end do
    
    !> @note Legacy commented code: Final mixing ratios for reaction amounts
    !> @details Would output final state of reaction amount mixing ratios at t=t_final
    !> Currently only outputs first element of each column
    ! write(unit,"(/,2x,'Final mixing ratios reaction amounts (by rows):'/)")
    ! do i=1,this%mixing_ratios_R%num_cols !< Loop over all columns for final reaction mixing ratios
    !     write(unit,"(2x,(ES15.5))") this%mixing_ratios_R%cols(i)%col_1(1) !< Output: only first R_k mixing ratio
    !     !> @note Nested commented code: Alternative full output
    !     !write(unit,"(2x,I5,*(ES15.5))") this%mixing_ratios_R%cols(i)%dim, (this%mixing_ratios_R%cols(i)%col_1(j), &
    !     !    j=1,this%mixing_ratios_R%cols(i)%dim) !< Would output: dim, then all final R_k mixing ratios
    ! end do
    
    !> @note Legacy commented code: Mixing waters indices
    !> @details Would output indices of water types contributing to mixing at each cell,
    !> plus counts of upstream and downstream waters for connectivity information
    ! write(unit,"(/,2x,'Mixing waters indices + no. upstream waters + no. downstream waters:'/)")
    ! do i=1,this%mix_conc_indices%num_cols !< Loop over all cells with mixing water information
    !     write(unit,"(2x,*(I5))") (this%mix_conc_indices%cols(i)%col_1(j), &
    !         j=1,this%mix_conc_indices%cols(i)%dim) !< Output: water indices, n_upstream, n_downstream [-]
    ! end do
    
    !> @note Legacy commented code: Reaction mixing waters indices
    !> @details Would output indices specifically for reactive mixing (dominant water types for reactions),
    !> plus connectivity information for reaction network
    ! write(unit,"(/,2x,'Reaction mixing waters indices + no. upstream waters + no. downstream waters:'/)")
    ! do i=1,this%mix_react_indices%num_cols !< Loop over all cells with dominant reaction water info
    !     write(unit,"(2x,*(I5))") (this%mix_react_indices%cols(i)%col_1(j), &
    !         j=1,this%mix_react_indices%cols(i)%dim) !< Output: dominant indices, n_upstream, n_downstream [-]
    ! end do
    
    !> @note Legacy commented code: Lagrangian particles information
    !> @details Would output Lagrangian particle data if Lagrangian tracking is enabled
    !> Lagrangian methods track water parcels as particles moving through the domain
    ! if (this%Lagr_flag) then !< Check if Lagrangian particle tracking is enabled (flag = .true.)
    !     !> @note Nested commented code: Lagrangian particle details
    !     !write(unit,"(/,2x,'Lagrangian particles information:'/)")  !< Section header
    !     !write(unit,"(2x,'Number of Lagrangian particles:',I5/)") this%num_Lagr_particles !< Total particle count
    !     !write(unit,"(2x,'Initial water particles coordinates (by rows):'/)")  !< Subsection header
    !     !do i=1,this%spatial_discr%num_targets !< Loop over all targets (particles initialized at target locations)
    !     !    write(unit,"(2x,*(ES15.5))") this%targets(i)%coord !< Output: particle x-coordinate [L]
    !     !end do
    ! end if
    close(unit)
    
end subroutine write_transport_data_WMA_2D !< End of write_transport_data_WMA subroutine