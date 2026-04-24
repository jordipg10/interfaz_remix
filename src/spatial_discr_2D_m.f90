!> \file spatial_discr_2D_m.f90
!> \brief Two-dimensional spatial discretization module.
!> \details
!> Defines homogeneous and heterogeneous 2D Eulerian mesh types extending
!> the 1D mesh classes. Adds a y-direction cell size \f$\Delta y\f$ and
!> computes cell diagonal lengths for stability analysis.
!>
!> \par Mesh Types:
!> - `mesh_2D_Euler_homog_c`   : uniform \f$\Delta x\f$ and \f$\Delta y\f$
!> - `mesh_2D_Euler_heterog_c` : variable cell sizes in both directions
!>
!> \see spatial_discr_m, spatial_discr_1D_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization

!> \brief 2D spatial discretization module.
module spatial_discr_2D_m
    use spatial_discr_m, only: spatial_discr_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
    save
    private
    !> \brief Homogeneous 2D Eulerian mesh with uniform cell sizes.
    !> \details Extends the 1D homogeneous mesh with a y-direction cell size.
    type, public, extends(mesh_1D_Euler_homog_c) :: mesh_2D_Euler_homog_c
        real(kind=8) :: Delta_y        !< [L] Uniform cell size in y direction
        real(kind=8) :: Delta_y_D      !< [-] Dimensionless cell size in y direction
        real(kind=8) :: sq_hypot       !< [L^2] Squared diagonal of cells \f$\Delta x^2 + \Delta y^2\f$
        integer(kind=4) :: Num_cells_x !< [#] Number of cells along x direction
        integer(kind=4) :: Num_cells_y !< [#] Number of cells along y direction
    contains
        procedure :: set_Delta_y_homog
        procedure :: read_mesh=>read_mesh_homog_2D
        procedure :: compute_Delta_x=>compute_Delta_x_2D_homog
        procedure :: get_Cell_size=>get_cell_size_2D_homog
        procedure :: get_max_cell_size=>get_max_cell_size_2D_homog
        procedure :: compute_measure=>compute_measure_homog_2D
        procedure :: compute_Delta_y
        procedure :: compute_sq_hypot_homog
        procedure :: compute_Num_targets=>compute_Num_targets_2D_homog
        procedure :: check_exit=>check_exit_2D_homog
        procedure :: refine_mesh=>refine_mesh_homog_2D
        procedure :: compute_dimless_mesh=>compute_dimless_mesh_2D_homog
        procedure :: compute_final_point=>compute_final_point_2D_homog
        procedure :: get_target_ind=>get_target_ind_2D_homog
        procedure :: get_num_cells=>get_num_cells_2D_homog
    end type
    
    !> \brief Heterogeneous 2D Eulerian mesh with variable cell sizes.
    !> \details Extends the 1D heterogeneous mesh with variable y-direction cell sizes.
    type, public, extends(mesh_1D_Euler_heterog_c) :: mesh_2D_Euler_heterog_c
        real(kind=8), allocatable :: Delta_y(:)     !< [L] Variable cell sizes in y direction
        real(kind=8), allocatable :: Delta_y_D(:)   !< [-] Dimensionless cell sizes in y direction
        real(kind=8), allocatable :: sq_hypot(:)    !< [L^2] Squared diagonals of cells
    contains
        procedure :: set_Delta_y_heterog
        procedure :: read_mesh=>read_mesh_heterog
        !procedure :: compute_Delta_x=>compute_Delta_x_2D_heterog
        procedure :: get_Cell_size=>get_cell_size_2D_heterog
        procedure :: get_max_cell_size=>get_max_cell_size_2D_heterog
        procedure :: compute_measure=>compute_measure_heterog_2D
        procedure :: check_exit=>check_exit_2D_heterog
        procedure :: refine_mesh=>refine_mesh_heterog_2D
        procedure :: compute_dimless_mesh=>compute_dimless_mesh_2D_heterog
        procedure :: get_target_ind=>get_target_ind_2D_heterog
        procedure :: compute_sq_hypot_heterog
        procedure :: compute_Num_targets=>compute_Num_targets_2D_heterog
        procedure :: get_num_cells=>get_num_cells_2D_heterog
    end type
    
    interface
        
        
        
    end interface
    
    contains
        subroutine set_Delta_y_homog(this,Delta_y)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            real(kind=8), intent(in) :: Delta_y
            this%Delta_y=Delta_y
        end subroutine
        
        subroutine set_Delta_y_heterog(this,Delta_y)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            real(kind=8), intent(in) :: Delta_y(:)
            this%Delta_y=Delta_y
        end subroutine
        
        subroutine read_mesh_homog_2D(this,filename,phi)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
            
            integer(kind=4) :: i, j !> Loop indices
            real(kind=8), allocatable :: coords(:) !> Coordinates of the targets (2D)
            character(len=100) :: line !> Line variable for reading file
            
            allocate(this%init_point(2),this%final_point(2),coords(2))
            open(unit=1,file=filename,status='old',action='read')
            do 
                read(1,*) line
                if (trim(line) == 'end') exit !> Exit loop if line reads 'end'
                if (trim(line) == 'SPATIAL DISCRETISATION') then
                    read(1,*) this%scheme
                    read(1,*) this%targets_flag
                    read(1,*) this%Delta_x
                    read(1,*) this%Delta_y
                    read(1,*) this%Num_cells_x
                    read(1,*) this%Num_cells_y
                    read(1,*) this%init_point
                    read(1,*) this%adapt_ref
                end if
            end do
            close(1)
            call this%set_dim(2) !> Set dimension to 2D
            this%Num_targets_defined=.true.
            !call this%compute_Delta_x()
            !call this%compute_Delta_y()
            call this%compute_Num_targets()
            call this%compute_measure()
            call this%allocate_targets()
            !> Boundary cells (including corners)
            do i=1,2
                coords(2)=this%init_point(2)+this%Delta_y*((2d0*((i-1)*this%Num_Cells_y+mod(i,2))-1d0)/2d0)
                do j=1,this%Num_cells_x
                    call this%targets((i-1)*this%Num_cells_x*(this%Num_cells_y-1)+j)%set_id(&
                        (i-1)*this%Num_cells_x*(this%Num_cells_y-1)+j)
                    call this%targets((i-1)*this%Num_cells_x*(this%Num_cells_y-1)+j)%set_boundary_flag(.true.)
                    call this%targets((i-1)*this%Num_cells_x*(this%Num_cells_y-1)+j)%set_measure(this%Delta_x*&
                        this%Delta_y*(1-this%targets_flag))
                    coords(1)=this%init_point(1)+this%Delta_x*((2d0*j-1d0)/2d0)
                    call this%targets((i-1)*this%Num_cells_x*(this%Num_cells_y-1)+j)%set_coordinates(coords)
                end do
            end do
            do j=1,2
                coords(1)=this%init_point(1)+this%Delta_x*((2d0*((j-1)*this%Num_cells_x+mod(j,2))-1d0)/2d0)
                do i=1,this%Num_cells_y-2
                    call this%targets((i+j-1)*this%Num_cells_x+mod(j,2))%set_id(&
                        (i+j-1)*this%Num_cells_x+mod(j,2))
                    call this%targets((i+j-1)*this%Num_cells_x+mod(j,2))%set_boundary_flag(.true.)
                    call this%targets((i+j-1)*this%Num_cells_x+mod(j,2))%set_measure(this%Delta_x*&
                        this%Delta_y*(1-this%targets_flag))
                    !coords(1)=this%init_point(1)+this%Delta_x*(1d0-1d0/(2d0-this%targets_flag))
                    coords(2)=this%init_point(2)+this%Delta_y*((2d0*(i+1)-1d0)/2d0)
                    call this%targets((i+j-1)*this%Num_cells_x+mod(j,2))%set_coordinates(coords)
                end do
            end do
            !> Non-boundary cells
            do i=2,this%Num_cells_y-1
                coords(2)=this%init_point(2)+this%Delta_y*((2d0*i-1d0)/2d0)
                do j=2,this%Num_cells_x-1
                    call this%targets((i-1)*this%Num_cells_x+j)%set_id((i-1)*this%Num_cells_x+j)
                    call this%targets((i-1)*this%Num_cells_x+j)%set_boundary_flag(.false.)
                    call this%targets((i-1)*this%Num_cells_x+j)%set_measure(&
                        this%Delta_x*this%Delta_y*(1-this%targets_flag))
                    coords(1)=this%init_point(1)+this%Delta_x*((2d0*j-1d0)/2d0)
                    call this%targets((i-1)*this%Num_cells_x+j)%set_coordinates(coords)
                end do
            end do
            ! !> Last target
            ! call this%targets(this%Num_targets)%set_id(this%Num_targets)
            ! call this%targets(this%Num_targets)%set_boundary_flag(.true.)
            ! call this%targets(this%Num_targets)%set_measure(&
            !     this%Delta_x*this%Delta_y*(1-this%targets_flag))
            ! coords(1)=this%init_point(1)+this%Delta_x*(this%Num_targets-1d0/(2d0-this%targets_flag))
            ! coords(2)=this%init_point(2)+this%Delta_y*(this%Num_targets-1d0/(2d0-this%targets_flag))
            ! call this%targets(this%Num_targets)%set_coordinates(coords)
            call this%compute_final_point()
        end subroutine
        
        subroutine read_mesh_heterog(this,filename,phi)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
            
            open(unit=1,file=filename,status='old',action='read')
            read(1,"(/,I10)") this%scheme
            read(1,*) this%targets_flag
            read(1,*) this%Num_targets
            allocate(this%Delta_x(this%Num_targets),this%Delta_y(this%Num_targets))
            read(1,*) this%Delta_x
            read(1,*) this%Delta_y
            close(1)
            call this%set_dim(2) !> Set dimension to 2D
            call this%compute_measure()
        end subroutine
        
        ! function get_Delta_y_homog(this,i) result(Delta_y)
        !     implicit none
        !     class(mesh_2D_Euler_homog_c) :: this
        !     integer(kind=4), intent(in), optional :: i
        !     real(kind=8) :: Delta_y
        !     Delta_y=this%Delta_y
        ! end function

        function get_cell_size_2D_homog(this,i) result(cell_size)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size
            cell_size=this%Delta_y*this%Delta_x
        end function
        
        ! function get_Delta_y_heterog(this,i) result(Delta_y)
        !     implicit none
        !     class(mesh_2D_Euler_heterog_c) :: this
        !     integer(kind=4), intent(in), optional :: i
        !     real(kind=8) :: Delta_y
        !     if (present(i)) then
        !         Delta_y=this%Delta_y(i) ! return cell size at index i
        !     else
        !         Delta_y=this%Delta_y(1) ! return minimum cell size if no index is provided
        !     end if
        ! end function

        function get_cell_size_2D_heterog(this,i) result(cell_size)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size
            if (present(i)) then
                cell_size=this%Delta_y(i)*this%Delta_x(i) ! return cell size at index i
            else
                cell_size=this%Delta_y(1)*this%Delta_x(1) ! return minimum cell size if no index is provided
            end if
        end function
        
        subroutine compute_measure_homog_2D(this)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%measure=(this%Num_targets-this%targets_flag)*this%Delta_x*this%Delta_y
        end subroutine
        
        subroutine compute_measure_heterog_2D(this)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            this%measure=sum(this%Delta_x*this%Delta_y)
        end subroutine
        
        subroutine compute_Delta_y(this)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%Delta_y=this%measure/(this%Delta_x*(this%Num_targets-this%targets_flag))
        end subroutine
        
        subroutine compute_Num_targets_2D_homog(this)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%Num_targets=(this%Num_cells_x + this%targets_flag)*(&
                this%Num_cells_y + this%targets_flag)
        end subroutine

        subroutine compute_Num_targets_2D_heterog(this)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            this%Num_targets=size(this%Delta_x) + this%targets_flag
        end subroutine
        
        function get_dim_2D_homog(this) result(dim)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            integer(kind=4) :: dim
            dim=2
        end function
        
        function get_dim_2D_heterog(this) result(dim)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            integer(kind=4) :: dim
            dim=2
        end function
        
        subroutine compute_dimless_mesh_2D_homog(this,char_length)
        implicit none
        class(mesh_2D_Euler_homog_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic measure for dimensionless form
        !if (this%dimless) then
            this%Delta_y_D=this%Delta_y/char_length
            this%measure_D=this%measure/char_length
            !this%init_point=this%init_point/char_measure
        !end if
        end subroutine
        
        subroutine compute_dimless_mesh_2D_heterog(this,char_length)
        implicit none
        class(mesh_2D_Euler_heterog_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic measure for dimensionless form
        this%Delta_y_D=this%Delta_y/char_length
        this%measure_D=this%measure/char_length
        !this%init_point=this%init_point/char_length
        end subroutine

        subroutine check_exit_2D_homog(this,coords,exit_flag)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit_flag !> exit flag
            
            exit_flag=.true.
            if (coords(1)<this%init_point(1) .or. coords(1)>this%final_point(1)) then
                print *, "Error: Coordinates out of bounds along x."
                !exit_flag=.true.
            else if (coords(2)<this%init_point(2) .or. coords(2)>this%final_point(2)) then
                print *, "Error: Coordinates out of bounds along y."
                !exit_flag=.true.
            else
                exit_flag=.false.
            end if
        end subroutine

        subroutine check_exit_2D_heterog(this,coords,exit_flag)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit_flag

            !exit_flag=.false.
            exit_flag=.true.
            if (coords(1)<this%init_point(1) .or. coords(1)>this%final_point(1)) then
                print *, "Error: Coordinates out of bounds along x."
                !exit_flag=.true.
            else if (coords(2)<this%init_point(2) .or. coords(2)>this%final_point(2)) then
                print *, "Error: Coordinates out of bounds along y."
            else
                exit_flag=.false.
            end if
        end subroutine

        function get_target_ind_2D_homog(this,coord) result(target_ind)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index
            
            integer(kind=4) :: i
            
            if (size(coord)/=2) error stop "Dimension error in get_target_ind_2D_homog"
            
            target_ind=0 ! Default value
            
            do i=1,this%Num_targets
                if ((abs(coord(1)-this%targets(i)%coord(1)) < this%Delta_x/2d0) .and. &
                    (abs(coord(2)-this%targets(i)%coord(2)) < this%Delta_y/2d0)) then
                    target_ind=i
                    return
                end if
            end do
            
        end function get_target_ind_2D_homog

        function get_target_ind_2D_heterog(this,coord) result(target_ind)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index

            integer(kind=4) :: target_ind_x !> target index
            integer(kind=4) :: target_ind_y !> target index
            integer(kind=4) :: i

            if (size(coord)/=2) error stop "Dimension error in get_target_ind_2D_heterog"

            target_ind=0 ! Default value

            if (coord(1) < this%init_point(1) .or. coord(1) > this%final_point(1)) then
                error stop "Error: Coordinates out of bounds along x."
            else if (coord(2) < this%init_point(2) .or. coord(2) > this%final_point(2)) then
                error stop "Error: Coordinates out of bounds along y."
            else
                !> Proceed to find target index along x direction
                if (coord(1) <= this%init_point(1) + 0.5*this%Delta_x(1)) then
                    target_ind_x = 1 ! Set target index for the first target
                end if
                if (coord(1) > this%final_point(1) - 0.5*this%Delta_x(this%Num_targets)) then
                    target_ind_x = this%Num_targets ! Set target index for the last target
                end if
                !> Proceed to find target index along y direction
                if (coord(2) <= this%init_point(2) + 0.5*this%Delta_y(1)) then
                    target_ind_y = 1 ! Set target index for the first target
                end if
                if (coord(2) > this%final_point(2) - 0.5*this%Delta_y(this%Num_targets)) then
                    target_ind_y = this%Num_targets ! Set target index for the last target
                end if
                !> Loop through targets to find correct index
                do i=2,this%Num_targets-1
                    if (coord(1)<=this%targets(i-1)%coord(1)+0.5*this%Delta_x(i-1)) then
                        target_ind_x = i-1 ! Set target index
                    else if (coord(1)<=this%targets(i)%coord(1)) then
                        target_ind_x = i ! Set target index
                    end if
                    if (coord(2)<=this%targets(i-1)%coord(2)+0.5*this%Delta_y(i-1)) then
                        target_ind_y = i-1 ! Set target index
                    else if (coord(2)<=this%targets(i)%coord(2)) then
                        target_ind_y = i ! Set target index
                    end if
                end do
                ! if (coord(1) <= this%final_point(1)) then
                !     target_ind = this%Num_targets ! Set target index for the last target
                ! end if
            end if
        end function get_target_ind_2D_heterog
        
        function get_max_cell_size_2D_heterog(this) result(max_cell_size)
        class(mesh_2D_Euler_heterog_c), intent(in) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=maxval(this%Delta_y*this%Delta_x)
        end function
        
        function get_max_cell_size_2D_homog(this) result(max_cell_size)
        class(mesh_2D_Euler_homog_c), intent(in) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=this%Delta_y*this%Delta_x
        end function

        subroutine refine_mesh_heterog_2D(this,conc,conc_ext,rel_tol)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
            integer(kind=4) :: j,n
            real(kind=8), allocatable :: Delta_y_new(:)
            
            n=this%Num_targets
            
            do j=1,n-1
                if (inf_norm_vec_real(conc(:,j)-conc(:,j+1))/inf_norm_vec_real(conc(:,j))>=rel_tol) then
                    Delta_y_new=[Delta_y_new,this%Delta_y(j)/2,this%Delta_y(j)/2]
                else
                    Delta_y_new=[Delta_y_new,this%Delta_y(j)]
                end if
            end do
            deallocate(this%Delta_y)
            this%Delta_y=Delta_y_new
            this%Num_targets=size(this%Delta_y)
        end subroutine

        subroutine refine_mesh_homog_2D(this,conc,conc_ext,rel_tol)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance

            real(kind=8), allocatable :: conc_ref(:,:)
            real(kind=8), allocatable :: conc_ext_ref(:,:)
    
            integer(kind=4) :: j,n,ratio
            real(kind=8) :: Delta_y_old,Delta_y_new
            
            n=this%Num_targets
            Delta_y_old=this%Delta_y
            Delta_y_new=Delta_y_old
            do j=1,n-1
                if (inf_norm_vec_real(conc(:,j)-conc(:,j+1))/inf_norm_vec_real(conc(:,j))>=rel_tol) then
                    Delta_y_new=Delta_y_old/2
                    Delta_y_old=Delta_y_new
                end if
            end do
            if (Delta_y_new<this%Delta_y) then
                ratio=nint(this%Delta_y/Delta_y_new)
                this%Num_targets=this%Num_targets*ratio
                this%Delta_y=Delta_y_new
                allocate(conc_ref(size(conc,1),this%Num_targets),conc_ext_ref(size(conc_ext,1),this%Num_targets))
                do j=1,this%Num_targets
                    conc_ref(:,j)=conc(:,ceiling(j/(1d0*ratio)))
                    conc_ext_ref(:,j)=conc_ext(:,ceiling(j/(1d0*ratio)))
                end do
                deallocate(conc,conc_ext)
                conc=conc_ref
                conc_ext=conc_ext_ref
            end if
        end subroutine

        subroutine compute_final_point_2D_homog(this)
            !> Compute final point based on initial point and cell sizes
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%final_point(1)=this%init_point(1)+&
                (this%Num_cells_x)*this%Delta_x
            this%final_point(2)=this%init_point(2)+&
                (this%Num_cells_y)*this%Delta_y
        end subroutine

        subroutine compute_sq_hypot_homog(this)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%sq_hypot=this%Delta_x**2+this%Delta_y**2
        end subroutine

        subroutine compute_sq_hypot_heterog(this)
            implicit none
            class(mesh_2D_Euler_heterog_c) :: this
            integer(kind=4) :: i
            allocate(this%sq_hypot(this%Num_targets-this%targets_flag))
            do i=1,this%Num_targets-this%targets_flag
                this%sq_hypot(i)=this%Delta_x(i)**2+this%Delta_y(i)**2
            end do
        end subroutine

        subroutine compute_Delta_x_2D_homog(this)
            implicit none
            class(mesh_2D_Euler_homog_c) :: this
            this%Delta_x=this%measure/(this%Delta_y*(this%Num_targets-this%targets_flag))
        end subroutine
        
        function get_num_cells_2D_homog(this,dim) result(num_cells) !< Get number of cells along a given dimension
            implicit none
            class(mesh_2D_Euler_homog_c), intent(in) :: this !> spatial discretisation object
            integer(kind=4), intent(in), optional :: dim !> spatial dimension (1, 2, or 3), optional for future use
            integer(kind=4) :: num_cells !> number of cells along the specified dimension
            if (present(dim)) then
                if (dim<1 .or. dim>3) then
                    error stop "Invalid dimension in get_num_cells"
                else
                    !select type (this)
                    !type is (mesh_1D_Euler_homog_c)
                    !        num_cells=this%Num_cells_x
                    !class is (mesh_2D_Euler_homog_c)
                        if (dim==1) then
                            num_cells=this%Num_cells_x
                        else if (dim==2) then
                            num_cells=this%Num_cells_y
                        end if
                    !class default
                    !    num_cells=this%Num_targets-this%targets_flag !> Default: return total number of cells (targets minus nodes if targets_flag=1)
                    !end select
                end if
            else
                num_cells=this%Num_targets-this%targets_flag !> If dim not specified, return total number of cells (targets minus nodes if targets_flag=1)
            end if
        end function get_num_cells_2D_homog
        
        function get_num_cells_2D_heterog(this,dim) result(num_cells) !< Get number of cells along a given dimension
            implicit none
            class(mesh_2D_Euler_heterog_c), intent(in) :: this !> spatial discretisation object
            integer(kind=4), intent(in), optional :: dim !> spatial dimension (1, 2, or 3), optional for future use
            integer(kind=4) :: num_cells !> number of cells along the specified dimension
            if (present(dim)) then
                if (dim<1 .or. dim>3) then
                    error stop "Invalid dimension in get_num_cells"
                else
                    !if (dim==1) then
                    !    num_cells=this%Num_cells_x
                    !else if (dim==2) then
                    !    num_cells=this%Num_cells_y
                    !end if
                end if
            else
                num_cells=this%Num_targets-this%targets_flag !> If dim not specified, return total number of cells (targets minus nodes if targets_flag=1)
            end if
        end function get_num_cells_2D_heterog

end module spatial_discr_2D_m