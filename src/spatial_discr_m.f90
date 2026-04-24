!> \file spatial_discr_m.f90
!> \brief Spatial discretization module for numerical mesh management
!> \details
!>   This module provides the abstract base class for spatial discretization
!>   of computational domains in 1D, 2D, and 3D.
!>
!>   Spatial discretization types:
!>   - Target-based: Either cells (volumes) or nodes (points)
!>   - Adaptive: Support for mesh refinement
!>   - Dimensionless: Normalized coordinates for dimensional analysis
!>
!>   Discretization schemes:
!>   - 1: Traditional Finite Differences (TFD) - second-order accurate
!>   - 2: Proposed Finite Differences (PFD) - second-order accurate (Petchamé-Guerrero et al., 2024)
!>   - 3: Upwind - first-order, stable for advection-dominated flows
!>
!>   Key features:
!>   - Abstract interface for multiple mesh types
!>   - Target tracking (cells or nodes)
!>   - Mesh refinement capabilities
!>   - Dimensionless coordinate system
!>   - Particle exit detection
!>
!>   The module defines deferred procedures that must be implemented
!>   by concrete subclasses (e.g., mesh_1D_c, mesh_rad_c).
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2025
module spatial_discr_m
    use target_m, only: target_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
    save
    private
    !> \brief Abstract spatial discretization superclass
    !> \details
    !>   Base class for all spatial discretization implementations.
    !>   
    !>   Target types:
    !>   - targets_flag=0: Cells (control volumes)
    !>   - targets_flag=1: Nodes (grid points)
    !>   
    !>   Coordinate system:
    !>   - Physical: measure (length/area/volume)
    !>   - Dimensionless: measure_D (normalized)
    !>   - Bounded: [init_point, final_point]
    !>   
    !>   Schemes:
    !>   - 1: Traditional FD (TFD) - O(Δx²) accuracy, potential instability
    !>   - 2: Proposed FD (PFD) - O(Δx²) accuracy, potential instability
    !>   - 3: Upwind - O(Δx) accuracy, stable for advection
    !>   
    !>   Adaptive refinement:
    !>   - adapt_ref=0: Fixed mesh
    !>   - adapt_ref=1: Adaptive mesh with error-based refinement
    !>   
    !>   Must be extended by concrete implementations that define
    !>   dimension-specific behavior (1D, 2D, 3D, radial, etc.).
    type, public, abstract :: spatial_discr_c !> spatial discretisation abstract superclass
        integer(kind=4) :: Num_targets      !> number of targets (cells or nodes)
        logical :: Num_targets_defined      !> TRUE if Num_targets defined, FALSE otherwise
        integer(kind=4) :: targets_flag     !> 0: cells (control volumes), 1: nodes (grid points)
        real(kind=8) :: measure             !> physical measure: length (1D), area (2D), volume (3D) [length units]
        real(kind=8) :: measure_D           !> dimensionless measure (normalized by characteristic length)
        real(kind=8), allocatable :: init_point(:)          !> initial/left boundary point coordinates
        real(kind=8), allocatable :: final_point(:)         !> final/right boundary point coordinates
        integer(kind=4) :: scheme           !> Spatial discretisation scheme: 1=CFD, 2=IFD, 3=Upwind
        integer(kind=4) :: adapt_ref        !> adaptive refinement flag: 0=NO, 1=YES
        type(target_c), allocatable :: targets(:) !> targets array (cells or nodes) of size Num_targets
        integer(kind=4) :: dim !> spatial dimension (1, 2, or 3)
    contains
        procedure :: set_targets_flag               !< Set whether targets are cells or nodes
        procedure :: set_dim                        !< Set spatial dimension
        procedure :: set_targets                    !< Set targets array
        procedure :: allocate_targets               !< Allocate memory for targets array
        procedure :: set_Num_targets                !< Set number of targets with validation
        procedure :: set_measure                    !< Set physical domain measure
        procedure :: set_scheme                     !< Set numerical scheme (CFD/IFD/Upwind)
        !procedure :: get_target_ind                !< [DISABLED] Get target index from coordinates
        procedure(read_mesh), public, deferred :: read_mesh !< Read mesh from file (dimension-specific)
        procedure(get_cell_size), public, deferred :: get_cell_size !< Get cell size (homogeneous or heterogeneous)
        procedure(get_max_cell_size), public, deferred :: get_max_cell_size !< Get maximum cell size in mesh
        procedure(compute_dimless_mesh), public, deferred :: compute_dimless_mesh !< Compute dimensionless coordinates
        procedure(get_target_ind), public, deferred :: get_target_ind !< Get target index from spatial coordinates
        procedure(compute_measure), public, deferred :: compute_measure !< Compute total domain measure
        procedure(compute_Num_targets), public, deferred :: compute_Num_targets !< [DISABLED] Compute number of targets
        procedure(refine_mesh), public, deferred :: refine_mesh !< Adaptive mesh refinement based on solution error
        procedure(check_exit), public, deferred :: check_exit !< Check if coordinates are outside domain
        procedure(get_num_cells), public, deferred :: get_num_cells !< Get spatial dimension
    end type
        
    !> \brief Abstract interface definitions for deferred procedures
    !> \details
    !>   These interfaces must be implemented by concrete subclasses.
    !>   Each interface defines the signature for dimension-specific operations.
    abstract interface
        function get_num_cells(this,dim) result(num_cells)
            import spatial_discr_c
            class(spatial_discr_c), intent(in) :: this
            integer(kind=4), intent(in), optional :: dim
            integer(kind=4) :: num_cells
        end function
        !> \brief Get maximum cell size in mesh
        !> \param[in] this Spatial discretization object
        !> \return max_cell_size Maximum cell dimension [length units]
        !> \details
        !>   Returns the largest cell size in the mesh.
        !>   Used for stability analysis (CFL condition) and error estimation.
        function get_max_cell_size(this) result(max_cell_size)
        import spatial_discr_c
        class(spatial_discr_c), intent(in) :: this
        real(kind=8) :: max_cell_size
        end function
    
    
        !> \brief Compute dimensionless mesh coordinates
        !> \param[in,out] this        Spatial discretization object
        !> \param[in]     char_length Characteristic length for normalization [length units]
        !> \details
        !>   Transforms physical coordinates to dimensionless form:
        !>   \f[
        !>     x_D = \frac{x}{L_c}
        !>   \f]
        !>   where \f$L_c\f$ is the characteristic length.
        !>   Updates measure_D and dimensionless target positions.
        subroutine compute_dimless_mesh(this,char_length)
        import spatial_discr_c
        class(spatial_discr_c) :: this
        real(kind=8), intent(in) :: char_length !> characteristic length
        end subroutine
        
        !> \brief Read mesh from input file
        !> \param[in,out] this     Spatial discretization object
        !> \param[in]     filename Mesh data file path
        !> \param[in]     phi      Optional porosity value
        !> \details
        !>   Reads mesh data from file containing:
        !>   - Number of targets/cells
        !>   - Cell sizes (homogeneous or heterogeneous)
        !>   - Boundary point coordinates
        !>   - Optional: porosity field
        subroutine read_mesh(this,filename,phi)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi
        end subroutine

        !> \brief Get cell size
        !> \param[in] this Spatial discretization object
        !> \param[in] i    Optional cell index (for heterogeneous meshes)
        !> \return cell_size Cell dimension [length units]
        !> \details
        !>   Returns cell size:
        !>   - Homogeneous mesh: Constant Δx for all cells
        !>   - Heterogeneous mesh: Δx(i) for cell i
        function get_cell_size(this,i) result(cell_size)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size
        end function

        !> \brief Check if coordinates are outside domain
        !> \param[in]  this      Spatial discretization object
        !> \param[in]  coords    Coordinates to check [length units]
        !> \param[out] exit_flag TRUE if outside domain, FALSE otherwise
        !> \details
        !>   Determines if given coordinates have exited the computational domain.
        !>   Used in Lagrangian particle tracking to detect particle loss.
        !>   
        !>   Exit conditions:
        !>   - 1D: x < x_min or x > x_max
        !>   - 2D: Outside bounding box or polygon
        !>   - 3D: Outside bounding volume
        subroutine check_exit(this,coords,exit_flag)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            real(kind=8), intent(in) :: coords(:) !> coordinates to check
            logical, intent(out) :: exit_flag !> TRUE if coords are outside the domain, FALSE otherwise
        end subroutine

        !> \brief Get spatial dimension
        !> \param[in] this Spatial discretization object
        !> \return dim Dimension (1, 2, or 3)
        !> \details
        !>   Returns the spatial dimension of the discretization.
        function get_dim(this) result(dim)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4) :: dim
        end function
        
        !> \brief Compute total domain measure
        !> \param[in,out] this Spatial discretization object
        !> \details
        !>   Computes the total domain measure:
        !>   - 1D: Total length \f$ L = \sum_i \Delta x_i \f$
        !>   - 2D: Total area \f$ A = \sum_i \Delta A_i \f$
        !>   - 3D: Total volume \f$ V = \sum_i \Delta V_i \f$
        !>   
        !>   Updates this%measure with the computed value.
        subroutine compute_measure(this)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
        end subroutine
        
        !> \brief Adaptive mesh refinement
        !> \param[in,out] this     Spatial discretization object
        !> \param[in,out] conc     Concentration array [Num_species × Num_targets]
        !> \param[in,out] conc_ext External concentration array
        !> \param[in]     rel_tol  Relative tolerance for refinement criterion
        !> \details
        !>   Refines mesh adaptively based on solution error estimates.
        !>   
        !>   Algorithm:
        !>   1. Compute local error indicator: \f$ \eta_i = \frac{|c_i^{n+1} - c_i^n|}{|c_i^{n+1}|} \f$
        !>   2. Flag cells where \f$ \eta_i > \text{rel_tol} \f$
        !>   3. Subdivide flagged cells
        !>   4. Interpolate solution to refined mesh
        !>   5. Update Num_targets and target array
        !>   
        !>   Increases resolution in regions with steep gradients or reactions.
        subroutine refine_mesh(this,conc,conc_ext,rel_tol)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
            !integer(kind=4), intent(out) :: n_new
        end subroutine

        !> \brief Get target index from spatial coordinates
        !> \param[in] this   Spatial discretization object
        !> \param[in] coord  Space coordinates [length units]
        !> \return target_ind Target index (1 to Num_targets), 0 if not found
        !> \details
        !>   Locates which target (cell or node) contains the given coordinates.
        !>   
        !>   Search methods:
        !>   - 1D: Binary search or linear scan
        !>   - 2D/3D: Spatial indexing or brute force
        !>   
        !>   Returns 0 if coordinates are outside domain.
        !>   Used for particle tracking and point source placement.
        function get_target_ind(this,coord) result(target_ind)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> space coordinates
            integer(kind=4) :: target_ind !> target index
        end function

        subroutine compute_Num_targets(this)
            import spatial_discr_c
            implicit none
            class(spatial_discr_c) :: this
        end subroutine
    end interface
    
    contains

        !> \brief Set spatial dimension
        !> \param[in,out] this Spatial discretization object
        !> \param[in]     dim  Spatial dimension (1, 2, or 3)
        !> \details
        !>   Sets the spatial dimension for the discretization.
        !>   Must be called before allocating coordinate arrays.
        subroutine set_dim(this,dim)
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in) :: dim
            this%dim=dim
        end subroutine

        !> \brief Set targets flag (cells vs nodes)
        !> \param[in,out] this Spatial discretization object
        !> \param[in]     flag 0 for cells (control volumes), 1 for nodes (grid points)
        !> \details
        !>   Specifies whether targets represent:
        !>   - 0: Cells (control volumes) - for finite volume methods
        !>   - 1: Nodes (grid points) - for finite difference/element methods
        !>   
        !>   Affects how fluxes, gradients, and measures are computed.
        subroutine set_targets_flag(this,flag)
            implicit none
            class(spatial_discr_c) :: this
            !integer(kind=4), intent(in) :: Num_targets
            integer(kind=4), intent(in) :: flag
            !this%Num_targets=Num_targets
            if (flag>1 .or. flag<0) error stop "Error in set_targets_flag"
            this%targets_flag=flag
        end subroutine
        
        !> \brief Set number of targets with validation
        !> \param[in,out] this        Spatial discretization object
        !> \param[in]     Num_targets Number of targets (must be positive)
        !> \details
        !>   Sets the number of computational targets (cells or nodes).
        !>   Validates that Num_targets > 0.
        !>   Must be called before allocating targets array.
        subroutine set_Num_targets(this,Num_targets)
            implicit none
            class(spatial_discr_c) :: this
            integer(kind=4), intent(in) :: Num_targets
            if (Num_targets<1) then
                error stop "Number of targets must be positive"
            end if
            this%Num_targets=Num_targets
        end subroutine 

        !> \brief Set targets array
        !> \param[in,out] this    Spatial discretization object
        !> \param[in]     targets Array of target objects
        !> \details
        !>   Assigns the targets array containing cell or node information.
        !>   Validates that array size matches Num_targets.
        !>   Each target contains: coordinates, volume/area, ID, etc.
        subroutine set_targets(this,targets)
            implicit none
            class(spatial_discr_c) :: this
            type(target_c), intent(in) :: targets(:)
            if (size(targets)/=this%Num_targets) then
                error stop "Dimension error in targets array"
            end if
            this%targets=targets
        end subroutine
        
        !> \brief Set physical domain measure
        !> \param[in,out] this    Spatial discretization object
        !> \param[in]     measure Physical domain measure [length, area, or volume units]
        !> \details
        !>   Sets the total domain measure:
        !>   - 1D: Total length L
        !>   - 2D: Total area A
        !>   - 3D: Total volume V
        !>   
        !>   Can be computed automatically via compute_measure() or set manually.
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
        
        !> \brief Set numerical discretization scheme
        !> \param[in,out] this   Spatial discretization object
        !> \param[in]     scheme Scheme number: 1=TFD, 2=PFD, 3=Upwind
        !> \details
        !>   Configures the spatial discretization scheme:
        !>   
        !>   1: Traditional Finite Differences (TFD)
        !>      - Second-order accurate: O(Δx²)
        !>      - Symmetric stencil
        !>      - May cause oscillations for advection-dominated flows
        !>      - Best for diffusion-dominated problems
        !>   
        !>   2: Proposed Finite Differences (PFD)
        !>      - Second-order accurate: O(Δx²)
        !>      - May cause oscillations for advection-dominated flows
        !>      - Best for diffusion-dominated problems
        !>   
        !>   3: Upwind
        !>      - First-order accurate: O(Δx)
        !>      - Stable for advection-dominated flows
        !>      - Asymmetric stencil (biased toward flow direction)
        !>      - Introduces numerical diffusion
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
        
        !> \brief Allocate memory for targets array
        !> \param[in,out] this Spatial discretization object
        !> \details
        !>   Allocates the targets array with size Num_targets.
        !>   If targets array already exists, it is deallocated first.
        !>   
        !>   Must call set_Num_targets() before this subroutine.
        !>   Each target will need to be initialized separately.
        subroutine allocate_targets(this)
        implicit none
        class(spatial_discr_c) :: this
        if (allocated(this%targets)) then
            deallocate(this%targets)
        end if
        allocate(this%targets(this%Num_targets))
        end subroutine allocate_targets
        
        !function get_num_cells(this,dim) result(num_cells) !< Get number of cells along a given dimension
        !    implicit none
        !    class(spatial_discr_c) :: this !> spatial discretisation object
        !    integer(kind=4), intent(in), optional :: dim !> spatial dimension (1, 2, or 3), optional for future use
        !    integer(kind=4) :: num_cells !> number of cells along the specified dimension
        !    if (present(dim)) then
        !        if (dim<1 .or. dim>3) then
        !            error stop "Invalid dimension in get_num_cells"
        !        else
        !            select type (this)
        !            !type is (mesh_1D_Euler_homog_c)
        !            !        num_cells=this%Num_cells_x
        !            class is (mesh_2D_Euler_homog_c)
        !                if (dim==1) then
        !                    num_cells=this%Num_cells_x
        !                else if (dim==2) then
        !                    num_cells=this%Num_cells_y
        !                end if
        !            !class default
        !            !    num_cells=this%Num_targets-this%targets_flag !> Default: return total number of cells (targets minus nodes if targets_flag=1)
        !            end select
        !        end if
        !    else
        !        num_cells=this%Num_targets-this%targets_flag !> If dim not specified, return total number of cells (targets minus nodes if targets_flag=1)
        !    end if
        !end function get_num_cells

end module
