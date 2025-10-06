!> Radial spatial discretisation module
module spatial_discr_rad_m
    use spatial_discr_m, only: spatial_discr_c
    implicit none
    save
    type, public, extends(spatial_discr_c) :: spatial_discr_rad_c
        !real(kind=8) :: thckn !> thickness of the radial mesh
        real(kind=8) :: r_max
        real(kind=8) :: r_max_D
        real(kind=8) :: r_min
        real(kind=8) :: r_min_D
        real(kind=8), allocatable :: Delta_r(:)
        real(kind=8), allocatable :: Delta_r_D(:)
    contains
        procedure, public :: compute_r_max
        procedure, public :: set_r_min
        procedure, public :: set_Delta_r
        procedure, public :: allocate_Delta_r
        procedure, public :: read_mesh=>read_mesh_rad_unif
        procedure, public :: get_Cell_size=>get_Delta_r
        procedure, public :: get_max_cell_size=>get_max_Delta_r
        !procedure, public :: get_dim=>get_dim_rad
        procedure, public :: set_r_max
        procedure, public :: compute_Delta_r
        procedure, public :: compute_measure=>compute_measure_rad
        procedure, public :: refine_mesh=>refine_mesh_rad
        procedure, public :: compute_dimless_mesh=>compute_dimless_mesh_rad
        procedure, public :: check_exit=>check_exit_rad
        procedure, public :: get_target_ind=>get_target_ind_rad
    end type
    
    interface
        
    end interface
    
    contains
        
        
        subroutine refine_mesh_rad(this,conc,conc_ext,rel_tol)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
        end subroutine
        
        subroutine set_r_min(this,r_min)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: r_min
        this%r_min=r_min
        end subroutine
        
        subroutine set_Delta_r(this,Delta_r)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: Delta_r(:)
            this%Delta_r=Delta_r
        end subroutine
        
        
        
        subroutine read_mesh_rad_unif(this,filename,phi)
        !> Reads radial mesh from a file
        !> We assume cell measure is uniform
            implicit none
            class(spatial_discr_rad_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi

            real(kind=8), parameter :: pi=4d0*atan(1d0) !> Pi constant
            real(kind=8) :: r_prev
            real(kind=8), allocatable :: r_i(:),coords(:)
            integer(kind=4) :: i

            !phi=0.01 !> autentica chapuza ad hoc
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%Num_targets
            read(1,*) this%dim !> Dimension of the mesh (1D, 2D, 3D)
            read(1,*) this%r_min
            !read(1,*) this%r_max
            read(1,*) this%measure !> Measure of the mesh (area in 2D, volume in 3D)
            read(1,*) this%adapt_ref
            close(1)
            this%Num_targets_defined=.true.
            !call this%set_dim(1) !> Set dimension (1D by default)
            call this%allocate_targets() !> Allocate targets array based on Num_targets
            call this%allocate_Delta_r() !> Allocate Delta_r array
            allocate(r_i(1)) !> Allocate r_i array for coordinates (chapuza)
            allocate(coords(1)) !> Allocate coordinates array
            !call this%targets(1)%set_id(1) !> Set ID for the first target
            !call this%targets(1)%set_boundary_flag(.true.) !> First target is a boundary target
            r_prev=this%r_min !> Initialize previous radius
            !coords(1)=this%r_min+(1d0-1d0/(2-this%targets_flag))*this%Delta_r(1) !> Compute coordinates for the first target
            !call this%targets(1)%set_coordinates(coords) !> Set coordinates for the first target
           ! call this%targets(1)%set_measure(pi*((this%r_min+this%Delta_r(1))**2-this%r_min**2)) !> Area in 2D, Volume in 3D (falta la porosidad y el espesor)
            do i=1,this%Num_targets
                call this%targets(i)%set_thickness() !> Set default thickness
                call this%targets(i)%set_measure((1-this%targets_flag)*this%measure/this%Num_targets) !> Set measure for each target
                r_i(1)=sqrt(r_prev**2+this%targets(i)%measure/(pi*this%targets(i)%thickness*phi)) !> Compute current radius
                this%Delta_r(i)=r_i(1)-r_prev !> Set Delta_r for the current target
                coords(1)=r_prev+(1d0-1d0/(2-this%targets_flag))*this%Delta_r(i) !> Compute coordinates of the target
                !call this%targets(i)%set_measure(pi*(r_i(1)**2-r_prev**2)) !> Area in 2D, Volume in 3D (falta la porosidad y el espesor)
                call this%targets(i)%set_coordinates(coords) !> Set coordinates for the current target
                call this%targets(i)%set_id(i) !> Set target ID
                call this%targets(i)%set_boundary_flag(.false.) !> Current target is not a boundary target
                r_prev=r_i(1) !> Update previous radius
            end do
            call this%targets(1)%set_boundary_flag(.true.) !> First target is a boundary target
            call this%targets(this%Num_targets)%set_boundary_flag(.true.) !> Last target is a boundary target
            call this%set_r_max(r_i(1)) !> Set maximum radius
        end subroutine
        
        function get_Delta_r(this,i) result(cell_size)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size !> Cell size (Delta_r)
            if (present(i)) then
                cell_size=this%Delta_r(i) !> Return specific Delta_r
            else
                cell_size=minval(this%Delta_r) !> Default to smallest Delta_r
            end if
        end function
        
        function get_dim_rad(this) result(dim)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4) :: dim
            dim=this%dim
        end function
        
        ! subroutine compute_r_max(this)
        !     implicit none
        !     class(spatial_discr_rad_c) :: this
        !     this%r_max=this%r_min+sum(this%Delta_r)
        ! end subroutine
        
        subroutine compute_Delta_r(this) !> Computes Delta_r based on r_min, r_max and Num_targets
            !! We assume uniform Delta_r
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4) :: i
            if (.not. allocated(this%Delta_r)) then
                allocate(this%Delta_r(this%Num_targets-this%targets_flag))
            end if
            do i=1,this%Num_targets-this%targets_flag
                this%Delta_r(i)=(this%r_max-this%r_min)/(this%Num_targets-this%targets_flag)
            end do
            !this%Delta_r=(this%r_max-this%r_min)/(this%Num_targets-this%targets_flag)
        end subroutine
        
        subroutine compute_measure_rad(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            if (this%dim == 1) then
                this%measure=this%r_max-this%r_min !> Length in 1D
            else if (this%dim == 2) then
                this%measure=pi*(this%r_max**2-this%r_min**2) !> Area in 2D
            else if (this%dim == 3) then
                this%measure=(4d0/3d0)*pi*(this%r_max**3-this%r_min**3) !> Volume in 3D
            else
                error stop "Dimension not implemented yet"
            end if
        end subroutine
        
        subroutine compute_dimless_mesh_rad(this,char_length)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic length scale for dimensionless conversion
        real(kind=8) :: r_max_dimless,r_min_dimless
        integer(kind=4) :: i
        this%r_max_D=this%r_max/char_length
        this%r_min_D=this%r_min/char_length
        this%Delta_r_D=this%Delta_r/char_length
        do i=1,this%Num_targets
            call this%targets(i)%compute_dimless_coords(char_length) !> Convert coordinates to dimensionless
        end do
        end subroutine

        subroutine allocate_Delta_r(this)
        implicit none
        class(spatial_discr_rad_c) :: this
        if (.not. allocated(this%Delta_r)) then
            allocate(this%Delta_r(this%Num_targets-this%targets_flag))
        end if
        end subroutine allocate_Delta_r

        subroutine compute_r_max(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            this%r_max=this%r_min+sum(this%Delta_r)
        end subroutine

        subroutine set_r_max(this,r_max)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: r_max
            if (r_max <= this%r_min) then
                error stop "r_max must be greater than r_min"
            end if
            this%r_max=r_max
        end subroutine

        subroutine check_exit_rad(this,coords,exit) !> Checks if the radial coordinate is within bounds
            implicit none
            class(spatial_discr_rad_c) :: this !> Radial spatial discretisation object
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit !> Exit flag to indicate if the check was successful

            exit = .false. ! Initialize exit flag to false

            if (coords(1) < this%r_min .or. coords(1) > this%r_max) then
                exit = .true.
                print *, "Error: Coordinates out of bounds."
            end if
            ! ! Check if the mesh is valid
            ! if (this%r_max <= this%r_min) then
            !     exit = .true.
            !     print *, "Error: r_max must be greater than r_min."
            !     return
            ! end if
            
            ! ! Check if Delta_r is allocated and has valid values
            ! if (.not. allocated(this%Delta_r)) then
            !     exit = .true.
            !     print *, "Error: Delta_r is not allocated."
            !     return
            ! end if
            
            ! if (any(this%Delta_r <= 0)) then
            !     exit = .true.
            !     print *, "Error: Delta_r must be positive."
            !     return
            ! end if
            
            !exit = .false. ! No errors found, continue execution
        end subroutine check_exit_rad

        function get_target_ind_rad(this,coord) result(target_ind)
            implicit none
            class(spatial_discr_rad_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index

            integer(kind=4) :: i
            real(kind=8), parameter :: eps=1d-12 !> small value for floating point comparison

            target_ind = 0 ! Initialize to 0 to indicate not found
            if (coord(1) < this%r_min .or. coord(1) > this%r_max ) then
                print *, "Error: Coordinates out of bounds."
                return
            end if
            if (coord(1) <= this%r_min + 0.5*this%Delta_r(1)) then
                target_ind = 1 ! Set target index for the first target
                return
            end if
            do i=2,this%Num_targets-1
                if (coord(1)<=this%targets(i-1)%coord(1)+0.5*this%Delta_r(i-1)) then
                    target_ind = i-1 ! Set target index
                    return
                else if (coord(1)<=this%targets(i)%coord(1)) then
                    target_ind = i ! Set target index
                    return
                else
                    continue ! Continue to next target
                end if
            end do
            if (coord(1) <= this%r_max) then
                target_ind = this%Num_targets ! Set target index for the last target
                return
            end if
        end function get_target_ind_rad
        
        function get_max_Delta_r(this) result(max_cell_size)
        class(spatial_discr_rad_c) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=maxval(this%Delta_r)
        end function

end module spatial_discr_rad_m
