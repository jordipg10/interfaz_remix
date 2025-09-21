!> 1D spatial discretisation module
module spatial_discr_1D_m
    use spatial_discr_m, only: spatial_discr_c
    implicit none
    save
    type, public, extends(spatial_discr_c) :: mesh_1D_Euler_homog_c
        real(kind=8) :: Delta_x !> Cell size
        real(kind=8) :: Delta_x_D !> Dimensionless cell size
    contains
        procedure, public :: set_Delta_x_homog
        procedure, public :: read_mesh=>read_mesh_homog
        procedure, public :: get_Cell_size=>get_Delta_x_homog
        procedure, public :: get_max_cell_size=>get_max_Delta_x_homog
        procedure, public :: compute_measure=>compute_measure_homog
        procedure, public :: compute_Delta_x
        procedure, private :: compute_Num_targets
        procedure, public :: check_exit=>check_exit_1D_homog
        procedure, public :: refine_mesh=>refine_mesh_homog
        procedure, public :: compute_dimless_mesh=>compute_dimless_mesh_1D_homog
        !procedure, public :: get_dim_1D=>get_dim_1D_homog
        procedure, public :: get_target_ind=>get_target_ind_1D_homog
    end type
    
    type, public, extends(spatial_discr_c) :: mesh_1D_Euler_heterog_c
        real(kind=8), allocatable :: Delta_x(:) !> Cell sizes
        real(kind=8), allocatable :: Delta_x_D(:) !> Dimensionless cell sizes
    contains
        procedure, public :: set_Delta_x_heterog
        procedure, public :: read_mesh=>read_mesh_heterog
        procedure, public :: get_Cell_size=>get_Delta_x_heterog
        procedure, public :: get_max_cell_size=>get_max_Delta_x_heterog
        procedure, public :: compute_measure=>compute_measure_heterog
        procedure, public :: check_exit=>check_exit_1D_heterog
        procedure, public :: refine_mesh=>refine_mesh_heterog
        procedure, public :: compute_dimless_mesh=>compute_dimless_mesh_1D_heterog
        procedure, public :: get_target_ind=>get_target_ind_1D_heterog
    end type
    
    interface
        subroutine refine_mesh_homog(this,conc,conc_ext,rel_tol)
            import mesh_1D_Euler_homog_c
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
        end subroutine
        
        subroutine refine_mesh_heterog(this,conc,conc_ext,rel_tol)
            import mesh_1D_Euler_heterog_c
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
        end subroutine
    end interface
    
    contains
        subroutine set_Delta_x_homog(this,Delta_x)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            real(kind=8), intent(in) :: Delta_x
            this%Delta_x=Delta_x
        end subroutine
        
        subroutine set_Delta_x_heterog(this,Delta_x)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            real(kind=8), intent(in) :: Delta_x(:)
            this%Delta_x=Delta_x
        end subroutine
        
        subroutine read_mesh_homog(this,filename,phi)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
            
            integer(kind=4) :: i !> Loop index
            real(kind=8) :: coords(1) !> Coordinates of the targets (1D)
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,*) this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%measure
            read(1,*) this%Num_targets
            read(1,*) this%init_point
            read(1,*) this%adapt_ref
            close(1)
            call this%set_dim(1) !> Set dimension to 1D
            this%Num_targets_defined=.true.
            call this%compute_Delta_x()
            call this%allocate_targets()
            !> First target
            call this%targets(1)%set_id(1)
            call this%targets(1)%set_boundary_flag(.true.)
            call this%targets(1)%set_measure(this%Delta_x*(1-this%targets_flag))
            coords=this%init_point+this%Delta_x*(1d0-1d0/(2d0-this%targets_flag))
            call this%targets(1)%set_coordinates(coords)
            !> Non-boundary targets
            do i=2,this%Num_targets-1
                call this%targets(i)%set_id(i)
                call this%targets(i)%set_boundary_flag(.false.)
                call this%targets(i)%set_measure(this%Delta_x*(1-this%targets_flag))
                coords=this%init_point+this%Delta_x*(i-1d0/(2d0-this%targets_flag))
                call this%targets(i)%set_coordinates(coords)
            end do
            !> Last target
            call this%targets(this%Num_targets)%set_id(this%Num_targets)
            call this%targets(this%Num_targets)%set_boundary_flag(.true.)
            call this%targets(this%Num_targets)%set_measure(this%Delta_x*(1-this%targets_flag))
            coords=this%init_point+this%Delta_x*(this%Num_targets-1d0/(2d0-this%targets_flag))
            call this%targets(this%Num_targets)%set_coordinates(coords)
        end subroutine
        
        subroutine read_mesh_heterog(this,filename,phi)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,I10)") this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%Num_targets
            allocate(this%Delta_x(this%Num_targets))
            read(1,*) this%Delta_x
            close(1)
            call this%set_dim(1) !> Set dimension to 1D
            call this%compute_measure()
        end subroutine
        
        function get_Delta_x_homog(this,i) result(Delta_x)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: Delta_x
            Delta_x=this%Delta_x
        end function
        
        function get_Delta_x_heterog(this,i) result(Delta_x)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: Delta_x
            if (present(i)) then
                Delta_x=this%Delta_x(i) ! return cell size at index i
            else
                Delta_x=minval(this%Delta_x) ! return minimum cell size if no index is provided
            end if
        end function
        
        subroutine compute_measure_homog(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%measure=(this%Num_targets-this%targets_flag)*this%Delta_x
        end subroutine
        
        subroutine compute_measure_heterog(this)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            this%measure=sum(this%Delta_x)
        end subroutine
        
        subroutine compute_Delta_x(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%Delta_x=this%measure/(this%Num_targets-this%targets_flag)
        end subroutine
        
        subroutine compute_Num_targets(this)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            this%Num_targets=this%measure/this%Delta_x + this%targets_flag
        end subroutine
        
        function get_dim_1D_homog(this) result(dim)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            integer(kind=4) :: dim
            dim=1
        end function
        
        function get_dim_1D_heterog(this) result(dim)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            integer(kind=4) :: dim
            dim=1
        end function
        
        subroutine compute_dimless_mesh_1D_homog(this,char_length)
        implicit none
        class(mesh_1D_Euler_homog_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic measure for dimensionless form
        !if (this%dimless) then
            this%Delta_x_D=this%Delta_x/char_length
            this%measure_D=this%measure/char_length
            !this%init_point=this%init_point/char_measure
        !end if
        end subroutine
        
        subroutine compute_dimless_mesh_1D_heterog(this,char_length)
        implicit none
        class(mesh_1D_Euler_heterog_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic measure for dimensionless form
        this%Delta_x_D=this%Delta_x/char_length
        this%measure_D=this%measure/char_length
        !this%init_point=this%init_point/char_length
        end subroutine

        subroutine check_exit_1D_homog(this,coords,exit)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit
            
            exit=.false.
            if (coords(1)<this%init_point .or. coords(1)>this%init_point+this%measure) then
                print *, "Error: Coordinates out of bounds."
                exit=.true.
            end if
        end subroutine

        subroutine check_exit_1D_heterog(this,coords,exit)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit

            exit=.false.
            if (coords(1)<this%init_point .or. coords(1)>this%init_point+this%measure) then
                print *, "Error: Coordinates out of bounds."
                exit=.true.
            end if
        end subroutine

        function get_target_ind_1D_homog(this,coord) result(target_ind)
            implicit none
            class(mesh_1D_Euler_homog_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index
            
            integer(kind=4) :: i
            
            if (size(coord)/=1) error stop "Dimension error in get_target_ind_1D_homog"
            
            target_ind=0 ! Default value
            
            do i=1,this%Num_targets
                if (abs(coord(1)-this%targets(i)%coord(1)) < this%Delta_x/2d0) then
                    target_ind=i
                    return
                end if
            end do
            
        end function get_target_ind_1D_homog

        function get_target_ind_1D_heterog(this,coord) result(target_ind)
            implicit none
            class(mesh_1D_Euler_heterog_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index

            integer(kind=4) :: i

            if (size(coord)/=1) error stop "Dimension error in get_target_ind_1D_heterog"

            target_ind=0 ! Default value

            if (coord(1) < this%init_point .or. coord(1) > this%init_point+this%measure) then
                print *, "Error: Coordinates out of bounds."
                return
            end if
            if (coord(1) <= this%init_point + 0.5*this%Delta_x(1)) then
                target_ind = 1 ! Set target index for the first target
                return
            end if
            do i=2,this%Num_targets-1
                if (coord(1)<=this%targets(i-1)%coord(1)+0.5*this%Delta_x(i-1)) then
                    target_ind = i-1 ! Set target index
                    return
                else if (coord(1)<=this%targets(i)%coord(1)) then
                    target_ind = i ! Set target index
                    return
                else
                    continue ! Continue to next target
                end if
            end do
            if (coord(1) <= this%init_point + this%measure) then
                target_ind = this%Num_targets ! Set target index for the last target
                return
            end if

        end function get_target_ind_1D_heterog
        
        function get_max_Delta_x_heterog(this) result(max_cell_size)
        class(mesh_1D_Euler_heterog_c) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=maxval(this%Delta_x)
        end function
        
        function get_max_Delta_x_homog(this) result(max_cell_size)
        class(mesh_1D_Euler_homog_c) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=this%Delta_x
        end function

end module
