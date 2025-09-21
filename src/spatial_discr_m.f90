!> Spatial discretisation module
!> This module defines the spatial discretisation class and its methods.
module spatial_discr_m
    use target_m, only: target_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
    save
    type, public, abstract :: spatial_discr_c !> spatial discretisation abstract superclass
        integer(kind=4) :: Num_targets      !> number of targets
        logical :: Num_targets_defined      !> TRUE if Num_targets defined, FALSE otherwise
        integer(kind=4) :: targets_flag     !> 0: cells
                                            !> 1: nodes
        real(kind=8) :: measure             !> length in 1D, area in 2D, volume in 3D
        real(kind=8) :: measure_D           !> dimensionless length in 1D, area in 2D, volume in 3D
        real(kind=8) :: init_point          !> initial point
        !real(kind=8) :: final_point         !> final point
        integer(kind=4) :: scheme           !> Spatial discretisation scheme:
                                                !> 1: CFD
                                                !> 2: IFD
                                                !> 3: Upwind
        integer(kind=4) :: adapt_ref        !> adaptive refinement (0: NO, 1: YES)
        type(target_c), allocatable :: targets(:) !> targets array (cells or nodes)
        integer(kind=4) :: dim !> dimension
    contains
        procedure, public :: set_targets_flag
        procedure, public :: set_dim
        procedure, public :: set_targets
        procedure, public :: allocate_targets
        procedure, public :: set_Num_targets
        procedure, public :: set_measure
        procedure, public :: set_scheme
        !procedure, public :: get_target_ind
        procedure(read_mesh), public, deferred :: read_mesh
        procedure(get_Cell_size), public, deferred :: get_Cell_size
        procedure(get_max_cell_size), public, deferred :: get_max_cell_size
        procedure(compute_dimless_mesh), public, deferred :: compute_dimless_mesh
        procedure(get_target_ind), public, deferred :: get_target_ind
        procedure(compute_measure), public, deferred :: compute_measure
        !procedure(compute_Num_targets), public, deferred :: compute_Num_targets
        procedure(refine_mesh), public, deferred :: refine_mesh
        procedure(check_exit), public, deferred :: check_exit        
    end type
        
    abstract interface
    
        function get_max_cell_size(this) result(max_cell_size)
        import spatial_discr_c
        class(spatial_discr_c) :: this
        real(kind=8) :: max_cell_size
        end function
    
    
        subroutine compute_dimless_mesh(this,char_length)
        import spatial_discr_c
        class(spatial_discr_c) :: this
        real(kind=8), intent(in) :: char_length !> characteristic length
        end subroutine
        
        subroutine read_mesh(this,filename,phi)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
        end subroutine
        
        function get_Cell_size(this,i) result(cell_size)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size
        end function

        subroutine check_exit(this,coords,exit)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            real(kind=8), intent(in) :: coords(:) !> coordinates to check
            logical, intent(out) :: exit
        end subroutine

        !> get dimension of the spatial discretisation
        function get_dim(this) result(dim)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4) :: dim
        end function
        
        subroutine compute_measure(this)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
        end subroutine
        
        subroutine refine_mesh(this,conc,conc_ext,rel_tol)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
            !integer(kind=4), intent(out) :: n_new
        end subroutine

        function get_target_ind(this,coord) result(target_ind)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index
        end function
    end interface
    
    contains

        subroutine set_dim(this,dim)
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in) :: dim
            this%dim=dim
        end subroutine

        subroutine set_targets_flag(this,flag)
            implicit none
            class(spatial_discr_c) :: this
            !integer(kind=4), intent(in) :: Num_targets
            integer(kind=4), intent(in) :: flag
            !this%Num_targets=Num_targets
            if (flag>1 .or. flag<0) error stop "Error in set_targets_flag"
            this%targets_flag=flag
        end subroutine
        
        subroutine set_Num_targets(this,Num_targets)
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in) :: Num_targets
            if (Num_targets<1) then
                error stop "Number of targets must be positive"
            end if
            this%Num_targets=Num_targets
        end subroutine 

        subroutine set_targets(this,targets)
            implicit none
            class(spatial_discr_c) :: this
            type(target_c), intent(in) :: targets(:)
            if (size(targets)/=this%Num_targets) then
                error stop "Dimension error in targets array"
            end if
            this%targets=targets
        end subroutine
        
        !subroutine read_mesh_1D(this,filename)
        !>    implicit none
        !>    class(spatial_discr_c) :: this
        !>    character(len=*), intent(in) :: filename
        !>    
        !>    real(kind=8), allocatable :: Delta_x(:)
        !>    
        !>    open(unit=1,file=filename,status='old',action='read')
        !>    read(1,*) this%Num_targets
        !>    select type (this)
        !>    type is (mesh_1D_Lagr_heterog_c)
        !>        allocate(Delta_x(this%Num_targets)) !> size(Delta_x_vec)=Num_targets
        !>        read(1,*) Delta_x
        !>        allocate(this%Delta_x(this%Num_targets)) 
        !>        this%Delta_x=Delta_x
        !>    end select
        !>    close(1)
        !>    this%Num_targets_defined=.true.
        !end subroutine 
        !
        !subroutine set_mesh_1D(this,Delta_x,Num_targets)
        !>    implicit none
        !>    class(spatial_discr_c) :: this
        !>    real(kind=8), intent(in) :: Delta_x
        !>    integer(kind=4), intent(in), optional :: Num_targets
        !>    select type (this)
        !>    type is (mesh_1D_Lagr_homog_c)
        !>        this%Delta_x=Delta_x
        !>        if (this%Num_targets_defined.eqv..false. .and. present(Num_targets)) then
        !>            this%Num_targets=Num_targets
        !>            this%Num_targets_defined=.true.
        !>        else if (this%Num_targets_defined.eqv..true. .and. present(Num_targets)) then
        !>            error stop "Num_targets already defined"
        !>        else if (this%Num_targets.eqv..false. .and. (.not. present(Num_targets))) then
        !>            error stop "Num_targets missing"
        !>        end if
        !>    end select
        !end subroutine
        !
        !subroutine set_mesh_1D_Lagr_homog(this,Delta_x,Num_targets)
        !>    implicit none
        !>    class(mesh_1D_Lagr_homog_c) :: this
        !>    real(kind=8), intent(in) :: Delta_x
        !>    integer(kind=4), intent(in), optional :: Num_targets
        !>    this%Delta_x=Delta_x
        !>    if (this%Num_targets_defined.eqv..false. .and. present(Num_targets)) then
        !>        this%Num_targets=Num_targets
        !>        this%Num_targets_defined=.true.
        !>    else if (this%Num_targets_defined.eqv..true. .and. present(Num_targets)) then
        !>        error stop "Num_targets already defined"
        !>    else if (this%Num_targets.eqv..false. .and. (.not. present(Num_targets))) then
        !>        error stop "Num_targets missing"
        !>    end if
        !end subroutine
        !
        !subroutine set_mesh_1D_Lagr_heterog(this,Delta_x,Num_targets)
        !>    implicit none
        !>    class(mesh_1D_Lagr_heterog_c) :: this
        !>    real(kind=8), intent(in) :: Delta_x(:)
        !>    integer(kind=4), intent(in), optional :: Num_targets
        !>    this%Delta_x=Delta_x
        !>    if (this%Num_targets_defined.eqv..true. .and. size(Delta_x)/=this%Num_targets) then
        !>        error stop "Dimension error in heterogeneous mesh"
        !>    else if (this%Num_targets_defined.eqv..false.) then
        !>        this%Num_targets=size(Delta_x)
        !>        this%Num_targets_defined=.true.
        !>    end if
        !end subroutine 
        
        !subroutine set_var_Delta_x(this,file_Delta_x)
        !>    implicit none
        !>    class(mesh_transport_1D) :: this
        !>    character(len=*), intent(in) :: file_Delta_x
        !>    
        !>    integer(kind=4) :: i,j,Num_elements
        !>    real(kind=8), allocatable :: Delta_x_vec(:)
        !>    !integer(kind=4), allocatable :: n_vec(:)
        !>    !> We assume file_Delta_x contains element sizes
        !>    select type (this)
        !>    type is (heterog_mesh_transport_1D)
        !>        open(unit=2,file=file_Delta_x,status='old',action='read')
        !>        read(2,*) Num_elements
        !>        allocate(Delta_x_vec(Num_elements)) !> size(Delta_x_vec)=Num_elements
        !>        read(2,*) Delta_x_vec
        !>        allocate(this%Delta_x(Num_elements)) 
        !>        this%Delta_x=Delta_x_vec
        !>        !print *, this%Delta_x
        !>        close(2)
        !>    class default
        !>        error stop "Wrong subclass"
        !>    end select
        !end subroutine set_var_Delta_x
        
        subroutine set_measure(this,measure)
            implicit none
            class(spatial_discr_c) :: this
            real(kind=8), intent(in) :: measure
            this%measure=measure
            !else if (this%targets_flag==0) then !> cells
            !    this%measure=this%get_Cell_size()
            !else !> nodes
            !    this%measure=0d0
            !end if
            !select type (this)
            !type is (homog_mesh_transport_1D)
            !>    this%measure=this%Num_elements*this%Delta_x
            !type is (heterog_mesh_transport_1D)
            !>    this%measure=sum(this%Delta_x)
            !end select
        end subroutine set_measure
        
        subroutine set_scheme(this,scheme)
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in) :: scheme
            if (scheme>3 .or. scheme<1) then
                error stop "Scheme not implemented yet"
            !else if (scheme.eqv.2 .and. this%targets_flag.eqv.0) then
                !error stop "Targets must be interfaces with IFDS"
            else
                this%scheme=scheme
            end if
        end subroutine 
        
        !function get_Cell_size(this) result(Delta_x)
        !>    implicit none
        !>    class(spatial_discr_c) :: this
        !>    real(kind=8), allocatable :: Delta_x(:)
        !>    select type (this)
        !>    type is (mesh_1D_Lagr_homog_c)
        !>        allocate(Delta_x(1))
        !>        Delta_x=get_Cell_size_homog(this)
        !>    type is (mesh_1D_Lagr_heterog_c)
        !>        Delta_x=get_Cell_size_heterog(this)
        !>    end select
        !end function
        !
        !function get_Cell_size_homog(this) result(Delta_x)
        !>    implicit none
        !>    class(mesh_1D_Lagr_homog_c) :: this
        !>    real(kind=8) :: Delta_x
        !>    Delta_x=this%Delta_x
        !end function
        !
        !function get_Cell_size_heterog(this) result(Delta_x)
        !>    implicit none
        !>    class(mesh_1D_Lagr_heterog_c) :: this
        !>    real(kind=8), allocatable :: Delta_x(:)
        !>    Delta_x=this%Delta_x
        !end function

        ! function get_target_ind(this,coord) result(target_ind)
        !     implicit none
        !     class(spatial_discr_c) :: this !> spatial discretisation class
        !     real(kind=8), intent(in) :: coord(:) !> space coordinates
        !     integer(kind=4) :: target_ind !> target index

        !     integer(kind=4) :: i
        !     real(kind=8), parameter :: eps=1d-12 !> small value for floating point comparison

        !     target_ind = 0 ! Initialize to 0 to indicate not found
        !     do i=2,this%Num_targets-1
        !         if (inf_norm_vec_real(this%targets(i-1)%coord+0.5*this%Delta_x(i-1))-coord)<eps) then
        !             target_ind = i
        !             return
        !         end if
        !     end do
        ! end function get_target_ind
        
        subroutine allocate_targets(this)
        implicit none
        class(spatial_discr_c) :: this
        if (allocated(this%targets)) then
            deallocate(this%targets)
        end if
        allocate(this%targets(this%Num_targets))
        end subroutine allocate_targets

end module
