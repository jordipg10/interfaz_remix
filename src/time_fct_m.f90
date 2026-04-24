!> \file time_fct_m.f90
!> \brief Time-dependent function module for piecewise-defined time series.
!> \details
!> Defines types for representing time-varying quantities (e.g., boundary
!> condition values that change at specified time step indices). Supports
!> integer and real-valued time series with I/O routines.
!>
!> \see time_discr_m, BCs_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization
module time_fct_m
implicit none
save
!> \brief Base time function class.
!> \details Stores a set of time indices and corresponding time values.
type, public :: time_fct_c
    integer(kind=4) :: ntime = 0              !< [#] Number of time entries
    integer(kind=4), allocatable :: time_ind(:)  !< [-] Time step indices
    real(kind=8), allocatable :: time(:)         !< [T] Time values
contains
    procedure :: set_ntime          !< Set number of time entries
    procedure :: allocate_time      !< Allocate index and value arrays
    procedure :: set_time           !< Set time values
    procedure :: set_time_ind       !< Set time step indices
    procedure :: get_time           !< Get time value at given index
end type time_fct_c

!> \brief Integer-valued time function subclass.
type, public, extends(time_fct_c) :: time_fct_int_c
    integer(kind=4), allocatable :: time_series(:)  !< [-] Integer time series values
contains
    procedure :: set_time_series=>set_time_series_int  !< Set integer time series
end type time_fct_int_c

!> \brief Real-valued time function subclass.
type, public, extends(time_fct_c) :: time_fct_real_c
    real(kind=8), allocatable :: time_series(:)  !< [-] Real time series values
contains
    procedure :: set_time_series=>set_time_series_real       !< Set real time series
    procedure :: read_time_series=>read_time_series_real     !< Read time series from file
    procedure :: allocate_time_series=>allocate_time_series_real  !< Allocate time series array
end type time_fct_real_c

    contains
    
    subroutine set_ntime(this,ntime)
        class(time_fct_c) :: this
        integer(kind=4), intent(in) :: ntime
        if (ntime < 0) then
            error stop "Error: ntime must be non-negative."
        end if
        this%ntime = ntime
    end subroutine 
    
    subroutine set_time(this,time)
        class(time_fct_c) :: this
        real(kind=8), intent(in) :: time(:)
        if (minval(time) < 0d0) then
            error stop "Error: attribute time must be non-negative."
        end if
        this%ntime = size(time)
        this%time = time
    end subroutine set_time
    
    function get_time(this,ind) result(time)
        class(time_fct_c), intent(in) :: this
        integer(kind=4), intent(in), optional :: ind !> Index of time step
        real(kind=8) :: time
    
        if (present(ind)) then
            if (ind < 1 .or. ind > this%ntime) then
                error stop "Error: Index out of bounds."
            end if
            time = this%time(ind)
        else
            time = this%time(1) ! Default to first time step
        end if
    end function get_time
    
    subroutine set_time_series_int(this,time_series)
        class(time_fct_int_c) :: this
        integer(kind=4), intent(in) :: time_series(:)
        this%ntime = size(time_series)
        this%time_series = time_series
        end subroutine set_time_series_int
    
        subroutine set_time_series_real(this,time_series)
        class(time_fct_real_c) :: this
        real(kind=8), intent(in) :: time_series(:)
        this%ntime = size(time_series)
        this%time_series = time_series
    end subroutine set_time_series_real
    
    subroutine read_time_series_real(this,unit)
    !> This subroutine reads a time series from a file and allocates the time series array.
    !> It assumes the first line of the file contains the number of time steps, 
    !! followed by the time series values.
    !> The file should be opened before calling this subroutine.
        class(time_fct_real_c) :: this !> Instance of time_fct_real_c
        integer(kind=4), intent(in) :: unit !> File unit number
        !character(len=*), intent(in) :: filename
    
        integer(kind=4) :: i !> Loop index
        !open(unit=1, file=filename, status='old', action='read') !> Open the file for reading
        read(unit,*) this%ntime !> Read number of time steps
        call this%allocate_time() !> Allocate time arrays
        call this%allocate_time_series() !> Allocate time series array
        do i = 1, this%ntime !> Read each time step
            read(unit,*) this%time_ind(i), this%time_series(i) !> Read time step index and value
        end do
        !close(1) !> Close the file
    end subroutine read_time_series_real
    
    subroutine allocate_time(this,ntime)
        class(time_fct_c) :: this
        integer(kind=4), intent(in), optional :: ntime
        if (present(ntime)) then
            if (ntime < 0) then
                error stop "Error: ntime must be non-negative."
            end if
            this%ntime = ntime
        end if
        allocate(this%time(this%ntime),this%time_ind(this%ntime))
    end subroutine allocate_time
    
    subroutine allocate_time_series_real(this,ntime)
        class(time_fct_real_c) :: this
        integer(kind=4), intent(in), optional :: ntime
        if (present(ntime)) then
            if (ntime < 0) then
                error stop "Error: ntime must be non-negative."
            end if
            this%ntime = ntime
        end if
        allocate(this%time_series(this%ntime))
    end subroutine allocate_time_series_real
    
    subroutine set_time_ind(this,ind)
    class(time_fct_c) :: this
    integer(kind=4), intent(in) :: ind(:)
    if (minval(ind) <= 0) then
        error stop "Error: Indices must be positive."
    else if (maxval(ind) > this%ntime) then
        error stop "Error: Indices exceed number of time steps."
    else
        this%time_ind = ind
    end if
    !allocate(ind(this%ntime))
    !ind = [(i, i = 1, this%ntime)]
    end subroutine set_time_ind
    
end module time_fct_m