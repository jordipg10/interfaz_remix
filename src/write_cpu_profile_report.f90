!> \file write_cpu_profile_report.f90
!> \brief Writes the CPU time profiling report to a file unit
!> \details
!>   Writes a formatted CPU profiling report containing:
!>   - Absolute elapsed times per phase (seconds and minutes)
!>   - Percentage breakdown relative to total subroutine time
!>   - Per-time-step and per-target-water averages for mixing and solver phases
!>
!>   Inner-loop phase timings are sampled (not measured every step).
!>   If sampling was used, the totals are scaled up to estimate full-run values
!>   before writing.
!>
!> \param[in] iunit_cpu            Fortran I/O unit number (file already open for writing)
!> \param[in] num_profiled_steps   Number of time steps where inner-loop profiling was active
!> \param[in] num_time_loc         Total number of time steps in the simulation
!> \param[in] num_non_can_vec      Number of non-canonical target waters (spatial points)
!> \param[in] cpu_total_elapsed    Total subroutine elapsed CPU time [s]
!> \param[in] cpu_update_old_total Accumulated CPU time for update_old_attributes [s]
!> \param[in] cpu_displace_total   Accumulated CPU time for particle displacement [s]
!> \param[in] cpu_data_prep_total  Accumulated CPU time for data preparation [s]
!> \param[in] cpu_mix_total        Accumulated CPU time for conservative mixing [s]
!> \param[in] cpu_getters_total    Accumulated CPU time for getter calls [s]
!> \param[in] cpu_solver_total     Accumulated CPU time for reactive mixing solver [s]
!> \param[in] cpu_derived_total    Accumulated CPU time for derived quantities [s]
!> \param[in] cpu_eq_rates_total   Accumulated CPU time for equilibrium reaction rates [s]
!> \param[in] cpu_solid_gas_total  Accumulated CPU time for solid/gas chemistry [s]
!> \param[in] cpu_output_total     Accumulated CPU time for output writing [s]
!> \param[in] cpu_snapshot_total   Accumulated CPU time for write_chemistry snapshots [s]
!> \param[in] dir                  Directory path for output files
!> \param[in] root                 Root name for output files
!> \author Jordi Petchamé-Guerrero
!> \date April 2026
!> \ingroup profiling
subroutine write_cpu_profile_report( &
        iunit_cpu, &
        num_profiled_steps, num_time_loc, num_non_can_vec, &
        cpu_total_elapsed, &
        cpu_update_old_total, cpu_displace_total, &
        cpu_data_prep_total, cpu_mix_total, cpu_getters_total, &
        cpu_solver_total, cpu_derived_total, cpu_eq_rates_total, &
        cpu_solid_gas_total, cpu_output_total, cpu_snapshot_total, &
        dir, root)
    implicit none
!> Arguments
    integer(kind=4), intent(in) :: iunit_cpu            !< Fortran I/O unit for the profiling file
    integer(kind=4), intent(in) :: num_profiled_steps   !< Number of inner-loop profiled time steps
    integer(kind=4), intent(in) :: num_time_loc         !< Total number of time steps
    integer(kind=4), intent(in) :: num_non_can_vec      !< Number of non-canonical target waters
    real(kind=8), intent(in) :: cpu_total_elapsed       !< Total subroutine elapsed CPU time [s]
    real(kind=8), intent(in) :: cpu_update_old_total    !< CPU time for update_old_attributes [s]
    real(kind=8), intent(in) :: cpu_displace_total      !< CPU time for particle displacement [s]
    real(kind=8), intent(in) :: cpu_data_prep_total     !< CPU time for data preparation [s]
    real(kind=8), intent(in) :: cpu_mix_total           !< CPU time for conservative mixing [s]
    real(kind=8), intent(in) :: cpu_getters_total       !< CPU time for getter calls [s]
    real(kind=8), intent(in) :: cpu_solver_total        !< CPU time for reactive mixing solver [s]
    real(kind=8), intent(in) :: cpu_derived_total       !< CPU time for derived quantities [s]
    real(kind=8), intent(in) :: cpu_eq_rates_total      !< CPU time for equilibrium reaction rates [s]
    real(kind=8), intent(in) :: cpu_solid_gas_total     !< CPU time for solid/gas chemistry [s]
    real(kind=8), intent(in) :: cpu_output_total        !< CPU time for output writing [s]
    real(kind=8), intent(in) :: cpu_snapshot_total      !< CPU time for write_chemistry snapshots [s]
    character(len=*), intent(in) :: dir                 !< Directory path for output files
    character(len=*), intent(in) :: root                !< Root name for output files
!> Local variables
    real(kind=8) :: profile_scale           !< Scaling factor for sampled profiling totals
    real(kind=8) :: s_data_prep             !< Scaled data preparation time [s]
    real(kind=8) :: s_mix                   !< Scaled conservative mixing time [s]
    real(kind=8) :: s_getters               !< Scaled getter calls time [s]
    real(kind=8) :: s_solver                !< Scaled reactive mixing solver time [s]
    real(kind=8) :: s_derived               !< Scaled derived quantities time [s]
    real(kind=8) :: s_eq_rates              !< Scaled equilibrium rates time [s]
    real(kind=8) :: s_solid_gas             !< Scaled solid/gas chemistry time [s]
    real(kind=8) :: s_output                !< Scaled output writing time [s]

    !> Copy input totals into local scaled copies (non-sampled phases stay unchanged)
    s_data_prep = cpu_data_prep_total
    s_mix       = cpu_mix_total
    s_getters   = cpu_getters_total
    s_solver    = cpu_solver_total
    s_derived   = cpu_derived_total
    s_eq_rates  = cpu_eq_rates_total
    s_solid_gas = cpu_solid_gas_total
    s_output    = cpu_output_total

    !> Scale sampled inner-loop profiling totals to estimate full-run values
    if (num_profiled_steps > 0 .and. num_profiled_steps < num_time_loc) then
        !> Compute scaling factor: ratio of total steps to profiled steps
        profile_scale = dble(num_time_loc) / dble(num_profiled_steps)
        !> Scale each sampled phase to full-run estimate
        s_data_prep = s_data_prep * profile_scale
        s_mix       = s_mix       * profile_scale
        s_getters   = s_getters   * profile_scale
        s_solver    = s_solver    * profile_scale
        s_derived   = s_derived   * profile_scale
        s_eq_rates  = s_eq_rates  * profile_scale
        s_solid_gas = s_solid_gas * profile_scale
        s_output    = s_output    * profile_scale
    end if

    !> Write profiling report header
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU TIME PROFILING: solve_RT_ideal_lump_Lagr_stat_flux_2D'
    write(iunit_cpu,*) '  (Inner-loop phases sampled from ', num_profiled_steps, ' of ', num_time_loc, ' time steps)'
    write(iunit_cpu,*) '====================================================================================='

    !> Write absolute elapsed times per phase in seconds and minutes
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Update old attributes:  ', cpu_update_old_total, ' s  (', cpu_update_old_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Particle displacement:  ', cpu_displace_total,   ' s  (', cpu_displace_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Data preparation:       ', s_data_prep,          ' s  (', s_data_prep/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Conservative mixing:    ', s_mix,                ' s  (', s_mix/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Getters (conc_nc+c1_oo):', s_getters,            ' s  (', s_getters/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Reactive mixing:        ', s_solver,             ' s  (', s_solver/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Derived quantities:     ', s_derived,            ' s  (', s_derived/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Equilibrium rates:      ', s_eq_rates,           ' s  (', s_eq_rates/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Solid/gas chemistry:    ', s_solid_gas,          ' s  (', s_solid_gas/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Output writing:         ', s_output,             ' s  (', s_output/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Snapshot (write_chem):  ', cpu_snapshot_total,   ' s  (', cpu_snapshot_total/60d0, ' min)'

    !> Write total subroutine time
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed,    ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='

    !> Write percentage breakdown if total elapsed time is positive
    if (cpu_total_elapsed > 0.0d0) then
        write(iunit_cpu,*) ''
        write(iunit_cpu,*) '  Breakdown (% of total):'
        write(iunit_cpu,'(A,F8.2,A)') '    Update old attributes:   ', 100d0*cpu_update_old_total/cpu_total_elapsed, ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Particle displacement:   ', 100d0*cpu_displace_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Data preparation:        ', 100d0*s_data_prep/cpu_total_elapsed,          ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Conservative mixing:     ', 100d0*s_mix/cpu_total_elapsed,                ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Getters (conc_nc+c1_oo): ', 100d0*s_getters/cpu_total_elapsed,            ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Reactive mixing:         ', 100d0*s_solver/cpu_total_elapsed,             ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Derived quantities:      ', 100d0*s_derived/cpu_total_elapsed,            ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Equilibrium rates:       ', 100d0*s_eq_rates/cpu_total_elapsed,           ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Solid/gas chemistry:     ', 100d0*s_solid_gas/cpu_total_elapsed,          ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Output writing:          ', 100d0*s_output/cpu_total_elapsed,             ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Snapshot (write_chem):   ', 100d0*cpu_snapshot_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,*) '====================================================================================='
    end if

    !> Write per-step averages
    write(iunit_cpu,*) ''
    write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing per time step:              ', s_mix / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing per time step/target water: ', s_mix / dble(num_time_loc * num_non_can_vec), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive mixing per time step:                  ', s_solver / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive mixing per time step/target water:     ', s_solver / dble(num_time_loc * num_non_can_vec), ' s'

    !> Close the CPU profiling output file
    close(iunit_cpu)
    !> Print notification that CPU profiling report has been written
    print *, 'CPU profiling written to: ', dir//root//'_cpu_profile.out'
end subroutine
