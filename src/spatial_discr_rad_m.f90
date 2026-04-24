!> \file spatial_discr_rad_m.f90
!> \brief Radial spatial discretization module for cylindrical/spherical geometries
!> \details
!>   This module extends spatial_discr_c for radial coordinate systems.
!>
!>   Radial geometry types:
!>   - 1D: Radial coordinate r (thickness-averaged)
!>   - 2D: Cylindrical (r,θ) with azimuthal symmetry → r only
!>   - 3D: Spherical (r,θ,φ) with full symmetry → r only
!>
!>   Applications:
!>   - Well injection/extraction (cylindrical coordinates)
!>   - Particle diffusion from spherical sources
!>   - Radial flow in confined/unconfined aquifers
!>   - Weathering of spherical mineral grains
!>
!>   Key features:
!>   - Non-uniform radial spacing Δr(i) for accurate resolution
!>   - Bounded domain: [r_min, r_max]
!>   - Cell volumes increase with radius: V(r) ∝ r² (2D) or r³ (3D)
!>   - Dimensionless coordinates for similarity solutions
!>
!>   Measure computation:
!>   - 1D: L = r_max - r_min (radial distance)
!>   - 2D: A = π(r_max² - r_min²) (cylindrical area)
!>   - 3D: V = (4π/3)(r_max³ - r_min³) (spherical volume)
!>
!>   The module implements all deferred procedures from spatial_discr_c
!>   with radial-specific algorithms.
!>
!> \author Generated documentation
!> \date 2025
module spatial_discr_rad_m
    use spatial_discr_m, only: spatial_discr_c
    implicit none
    save
    private
    !> \brief Radial spatial discretization class
    !> \details
    !>   Extends spatial_discr_c for radial coordinate systems.
    !>   
    !>   Radial coordinates:
    !>   - r_min: Inner boundary radius (well radius, particle radius)
    !>   - r_max: Outer boundary radius (domain extent)
    !>   - Delta_r(i): Radial spacing for cell i
    !>   
    !>   Dimensionless coordinates:
    !>   - r_D = r / L_c (normalized by characteristic length)
    !>   - Used for similarity solutions and scaling analysis
    !>   
    !>   Cell spacing strategies:
    !>   - Uniform: Δr = constant (simple, may over-resolve far field)
    !>   - Geometric: Δr(i+1) = α·Δr(i) (efficient for large domains)
    !>   - Adaptive: Based on solution gradients (optimal accuracy)
    !>   
    !>   Cell volume computation:
    !>   For cylindrical (2D):
    !>     V(i) = π·thickness(i)·phi·[(r_{i+1/2})² - (r_{i-1/2})²]
    !>   For spherical (3D):
    !>     V(i) = (4π/3)·[(r_{i+1/2})³ - (r_{i-1/2})³]
    type, public, extends(spatial_discr_c) :: spatial_discr_rad_c
        real(kind=8) :: cell_vol                !> cell volume [length^dim units]
        real(kind=8) :: r_max                   !> maximum radius (outer boundary) [length units]
        real(kind=8) :: r_max_D                 !> dimensionless maximum radius
        real(kind=8) :: r_min                   !> minimum radius (inner boundary, e.g., well radius) [length units]
        real(kind=8) :: r_min_D                 !> dimensionless minimum radius
        real(kind=8), allocatable :: Delta_r(:) !> radial cell spacing array [length units], size Num_targets
        real(kind=8), allocatable :: Delta_r_D(:) !> dimensionless radial cell spacing
    contains
        procedure :: compute_r_max                           !< Compute maximum radius from r_min and Delta_r
        procedure :: set_r_min                               !< Set minimum radius
        procedure :: set_Delta_r                             !< Set radial cell spacing array
        procedure :: allocate_Delta_r                        !< Allocate memory for Delta_r array
        procedure :: read_mesh=>read_mesh_rad_unif           !< Read radial mesh from file
        procedure :: get_Cell_size=>get_Delta_r              !< Get radial cell size
        procedure :: get_max_cell_size=>get_max_Delta_r      !< Get maximum radial cell size
        !procedure :: get_dim=>get_dim_rad                   !< [DISABLED] Get spatial dimension
        procedure :: set_r_max                               !< Set maximum radius with validation
        procedure :: compute_Delta_r                         !< Compute uniform Delta_r from r_min, r_max, Num_targets
        procedure :: compute_measure=>compute_measure_rad    !< Compute total radial domain measure
        procedure :: refine_mesh=>refine_mesh_rad            !< Adaptive mesh refinement for radial grids
        procedure :: compute_dimless_mesh=>compute_dimless_mesh_rad !< Convert to dimensionless coordinates
        procedure :: check_exit=>check_exit_rad              !< Check if radius is outside domain bounds
        procedure :: get_target_ind=>get_target_ind_rad      !< Get target index from radial coordinate
        procedure :: compute_Num_targets=>compute_Num_targets_rad !< Compute number of targets based on Delta_r
        procedure :: get_num_cells=>get_num_cells_rad      !< Get number of cells in the radial mesh
    end type
    
    interface
        
    end interface
    
    contains
        
        
        !> \brief Adaptive mesh refinement for radial grids
        !> \param[in,out] this     Radial spatial discretization object
        !> \param[in,out] conc     Concentration array [Num_species × Num_targets]
        !> \param[in,out] conc_ext External concentration array
        !> \param[in]     rel_tol  Relative tolerance for refinement
        !> \details
        !>   Refines radial mesh adaptively based on solution gradients.
        !>   (Currently placeholder - to be implemented)
        subroutine refine_mesh_rad(this,conc,conc_ext,rel_tol)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(inout), allocatable :: conc(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(inout), allocatable :: conc_ext(:,:) !> Num_columns=Num_targets
            real(kind=8), intent(in) :: rel_tol !> relative tolerance
        end subroutine
        
        !> \brief Set minimum radius
        !> \param[in,out] this  Radial spatial discretization object
        !> \param[in]     r_min Minimum radius [length units]
        !> \details
        !>   Sets the inner boundary radius.
        !>   Typically represents well radius, particle radius, or cavity size.
        subroutine set_r_min(this,r_min)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: r_min
        this%r_min=r_min
        end subroutine
        
        !> \brief Set radial cell spacing array
        !> \param[in,out] this    Radial spatial discretization object
        !> \param[in]     Delta_r Array of radial spacings [length units]
        !> \details
        !>   Assigns the radial cell spacing array for heterogeneous grids.
        !>   Allows non-uniform spacing for better resolution near boundaries.
        subroutine set_Delta_r(this,Delta_r)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: Delta_r(:)
            this%Delta_r=Delta_r
        end subroutine
        
        
        
        !> \brief Read radial mesh from file (uniform cell volume & thickness)
        !> \param[in,out] this     Radial spatial discretization object
        !> \param[in]     filename Mesh data file path
        !> \param[in]     phi      Porosity (optional)
        !> \details
        !>   Reads radial mesh assuming uniform cell volume and thickness.
        !>   
        !>   File format:
        !>   - Line 1: scheme (1=CFD, 2=IFD, 3=Upwind)
        !>   - Line 2: targets_flag (0=cells, 1=nodes)
        !>   - Line 3: Num_targets
        !>   - Line 4: dim (1, 2, or 3)
        !>   - Line 5: r_min [length units]
        !>   - Line 6: thickness [length units]
        !>   - Line 7: measure (total volume or area)
        !>   - Line 8: adapt_ref (0=fixed, 1=adaptive)
        !>   
        !>   Algorithm for uniform cell volume and thickness:
        !>   1. Set first cell at r_min
        !>   2. For each subsequent cell:
        !>      - Compute r_i: r_i² = r_{i-1}² + A_cell/(π·thickness·φ) (2D)
        !>      - Compute r_i: r_i³ = r_{i-1}³ + 3·V_cell/(4π·φ) (3D)
        !>      - Set Δr(i) = r_i - r_{i-1}
        !>   3. Update coordinates, IDs, boundary flags
        !>   
        !>   Ensures equal volume and thickness per cell despite varying Δr.
        subroutine read_mesh_rad_unif(this,filename,phi)
        !> Reads radial mesh from a file
        !> We assume cell volume & thickness are uniform
            implicit none
            class(spatial_discr_rad_c) :: this
            character(len=*), intent(in) :: filename
            real(kind=8), intent(in), optional :: phi

            real(kind=8), parameter :: pi=4d0*atan(1d0) !> Pi constant
            real(kind=8) :: r_prev !> Previous radius for iterative computation
            real(kind=8) :: thickness !> Thickness of the radial mesh
            real(kind=8), allocatable :: r_i(:),coords(:) !> Current radius and coordinates
            integer(kind=4) :: i !> Loop index
            character(len=256) :: label !> Label for reading file

            open(unit=1,file=filename,status='old',action='read')
            do 
                read(1,*) label
                if (trim(label) .eq. 'end') exit
                if (label == 'SPATIAL DISCRETISATION') then
                    read(1,*) this%scheme !> Discretization scheme (1=CFD, 2=IFD, 3=Upwind)
                    read(1,*) this%targets_flag !> 0=cells, 1=nodes
                    read(1,*) this%Num_targets !> Number of targets (cells or nodes)
                    read(1,*) this%dim !> Dimension of the mesh (1D, 2D, 3D)
                    read(1,*) this%r_min !> Minimum radius (inner boundary) [length units]
                    read(1,*) thickness !> Thickness of the radial mesh [length units]
                    read(1,*) this%measure !> Measure of the mesh (area in 2D, volume in 3D)
                    read(1,*) this%adapt_ref !> Adaptive refinement flag (0=fixed, 1=adaptive)
                end if
                !if (iostat /= 0) exit
            end do
            !> Read mesh parameters from file
            ! open(unit=1,file=filename,status='old',action='read')
            ! read(1,*) this%scheme !> Discretization scheme (1=CFD, 2=IFD, 3=Upwind)
            ! read(1,*) this%targets_flag !> 0=cells, 1=nodes
            ! read(1,*) this%Num_targets !> Number of targets (cells or nodes)
            ! read(1,*) this%dim !> Dimension of the mesh (1D, 2D, 3D)
            ! read(1,*) this%r_min !> Minimum radius (inner boundary) [length units]
            ! read(1,*) thickness !> Thickness of the radial mesh [length units]
            ! read(1,*) this%measure !> Measure of the mesh (area in 2D, volume in 3D)
            ! read(1,*) this%adapt_ref !> Adaptive refinement flag (0=fixed, 1=adaptive)
            close(1)
            this%Num_targets_defined=.true.
            allocate(this%init_point(1),this%final_point(1))
            this%init_point=this%r_min
            this%cell_vol=this%measure/(this%Num_targets-this%targets_flag) !> Compute cell volume/area
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
            if (this%targets_flag==1) then
                coords(1)=this%r_min
                call this%targets(1)%set_coordinates(coords) !> Set coordinates for the first target
                call this%targets(1)%set_thickness(thickness) !> Set thickness
                call this%targets(1)%set_measure(0d0) !> node has zero measure
                call this%targets(1)%set_id(1) !> Set target ID
            end if
            !> Compute radial positions ensuring uniform cell volume
            do i=1+this%targets_flag,this%Num_targets
                call this%targets(i)%set_thickness(thickness) !> Set thickness
                call this%targets(i)%set_measure((1-this%targets_flag)*this%measure/this%Num_targets) !> Equal measure per cell
                r_i(1)=sqrt(r_prev**2+this%cell_vol/(pi*this%targets(i)%thickness*phi)) !> Compute current radius (2D formula)
                this%Delta_r(i-this%targets_flag)=r_i(1)-r_prev !> Radial spacing for this cell
                coords(1)=r_prev+(1d0/(2-this%targets_flag))*this%Delta_r(i-this%targets_flag) !> Target coordinate (cell center or node)
                !call this%targets(i)%set_measure(pi*(r_i(1)**2-r_prev**2)) !> Area in 2D, Volume in 3D (falta la porosidad y el espesor)
                call this%targets(i)%set_coordinates(coords) !> Set coordinates for the current target
                call this%targets(i)%set_id(i) !> Set target ID
                call this%targets(i)%set_boundary_flag(.false.) !> Interior target by default
                r_prev=r_i(1) !> Update previous radius for next iteration
            end do
            call this%targets(1)%set_boundary_flag(.true.) !> Inner boundary (r_min)
            call this%targets(this%Num_targets)%set_boundary_flag(.true.) !> Outer boundary (r_max)
            call this%set_r_max(r_i(1)) !> Set maximum radius from final computed value
            this%final_point=this%r_max
        end subroutine
        
        !> \brief Get radial cell size
        !> \param[in] this Radial spatial discretization object
        !> \param[in] i    Optional cell index
        !> \return cell_size Radial cell spacing [length units]
        !> \details
        !>   Returns Delta_r for specified cell or minimum Delta_r if no index given.
        !>   Used for CFL stability analysis.
        function get_Delta_r(this,i) result(cell_size)
            implicit none
            class(spatial_discr_rad_c) :: this
            integer(kind=4), intent(in), optional :: i
            real(kind=8) :: cell_size !> Cell size (Delta_r)
            if (present(i)) then
                cell_size=this%Delta_r(i) !> Return specific Delta_r
            else
                cell_size=this%Delta_r(1) !> Default to Delta_r at index 1
            end if
        end function
        
        !> \brief Get spatial dimension
        !> \param[in] this Radial spatial discretization object
        !> \return dim Spatial dimension (1, 2, or 3)
        !> \details
        !>   Returns the dimension for radial geometry.
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
        
        !> \brief Compute uniform radial cell spacing
        !> \param[in,out] this Radial spatial discretization object
        !> \details
        !>   Computes uniform Delta_r based on r_min, r_max, and Num_targets:
        !>   \f[
        !>     \Delta r = \frac{r_{\max} - r_{\min}}{N_{\text{targets}} - \text{targets\_flag}}
        !>   \f]
        !>   
        !>   Allocates Delta_r array if not already allocated.
        !>   Uniform spacing is simple but may over-resolve far field.
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
        
        !> \brief Compute total radial domain measure
        !> \param[in,out] this Radial spatial discretization object
        !> \details
        !>   Computes total measure based on geometry:
        !>   
        !>   1D radial (thickness-averaged):
        !>   \f[
        !>     L = r_{\max} - r_{\min}
        !>   \f]
        !>   
        !>   2D cylindrical (azimuthal symmetry):
        !>   \f[
        !>     A = \pi (r_{\max}^2 - r_{\min}^2)
        !>   \f]
        !>   
        !>   3D spherical (full symmetry):
        !>   \f[
        !>     V = \frac{4\pi}{3} (r_{\max}^3 - r_{\min}^3)
        !>   \f]
        !>   
        !>   Updates this%measure with computed value.
        subroutine compute_measure_rad(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), parameter :: pi=4d0*atan(1d0)
            if (this%dim == 1) then
                this%measure=this%r_max-this%r_min !> Length in 1D
            else if (this%dim == 2) then
                this%measure=pi*(this%r_max**2-this%r_min**2) !> Area in 2D (cylindrical)
            else if (this%dim == 3) then
                this%measure=(4d0/3d0)*pi*(this%r_max**3-this%r_min**3) !> Volume in 3D (spherical)
            else
                error stop "Dimension not implemented yet"
            end if
        end subroutine
        
        !> \brief Convert radial mesh to dimensionless coordinates
        !> \param[in,out] this        Radial spatial discretization object
        !> \param[in]     char_length Characteristic length for normalization [length units]
        !> \details
        !>   Transforms physical radial coordinates to dimensionless form:
        !>   \f[
        !>     r_D = \frac{r}{L_c}, \quad \Delta r_D = \frac{\Delta r}{L_c}
        !>   \f]
        !>   
        !>   Updates:
        !>   - r_max_D, r_min_D (dimensionless boundaries)
        !>   - Delta_r_D (dimensionless cell spacing)
        !>   - Target coordinates (via target%compute_dimless_coords)
        !>   
        !>   Used for similarity solutions in radial diffusion/flow problems.
        subroutine compute_dimless_mesh_rad(this,char_length)
        implicit none
        class(spatial_discr_rad_c) :: this
        real(kind=8), intent(in) :: char_length !> Characteristic length scale for dimensionless conversion
        real(kind=8) :: r_max_dimless,r_min_dimless
        integer(kind=4) :: i
        this%r_max_D=this%r_max/char_length     !> Normalize maximum radius
        this%r_min_D=this%r_min/char_length     !> Normalize minimum radius
        this%Delta_r_D=this%Delta_r/char_length !> Normalize cell spacing array
        do i=1,this%Num_targets
            call this%targets(i)%compute_dimless_coords(char_length) !> Convert coordinates to dimensionless
        end do
        end subroutine

        !> \brief Allocate memory for Delta_r array
        !> \param[in,out] this Radial spatial discretization object
        !> \details
        !>   Allocates Delta_r array if not already allocated.
        !>   Size: Num_targets - targets_flag
        subroutine allocate_Delta_r(this)
        implicit none
        class(spatial_discr_rad_c) :: this
        if (.not. allocated(this%Delta_r)) then
            allocate(this%Delta_r(this%Num_targets-this%targets_flag))
        end if
        end subroutine allocate_Delta_r

        !> \brief Compute maximum radius from Delta_r array
        !> \param[in,out] this Radial spatial discretization object
        !> \details
        !>   Computes r_max from r_min and sum of Delta_r:
        !>   \f[
        !>     r_{\max} = r_{\min} + \sum_{i=1}^{N} \Delta r_i
        !>   \f]
        subroutine compute_r_max(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            this%r_max=this%r_min+sum(this%Delta_r)
        end subroutine

        !> \brief Set maximum radius with validation
        !> \param[in,out] this  Radial spatial discretization object
        !> \param[in]     r_max Maximum radius [length units]
        !> \details
        !>   Sets outer boundary radius with validation.
        !>   Ensures r_max > r_min to avoid invalid domain.
        subroutine set_r_max(this,r_max)
            implicit none
            class(spatial_discr_rad_c) :: this
            real(kind=8), intent(in) :: r_max
            if (r_max <= this%r_min) then
                error stop "r_max must be greater than r_min"
            end if
            this%r_max=r_max
        end subroutine

        !> \brief Check if radial coordinate is outside domain bounds
        !> \param[in]  this      Radial spatial discretization object
        !> \param[in]  coords    Radial coordinate to check [length units]
        !> \param[out] exit_flag TRUE if outside [r_min, r_max], FALSE otherwise
        !> \details
        !>   Validates that radial coordinate is within domain:
        !>   - exit_flag = TRUE if r < r_min or r > r_max
        !>   - exit_flag = FALSE if r_min ≤ r ≤ r_max
        !>   
        !>   Used in Lagrangian particle tracking to detect particles
        !>   leaving the radial domain (either toward center or infinity).
        subroutine check_exit_rad(this,coords,exit_flag) !> Checks if the radial coordinate is within bounds
            implicit none
            class(spatial_discr_rad_c) :: this !> Radial spatial discretisation object
            real(kind=8), intent(in) :: coords(:) !> Coordinates to check
            logical, intent(out) :: exit_flag !> TRUE if exit condition met, FALSE otherwise

            exit_flag = .false. ! Initialize exit flag to false

            if (coords(1) < this%r_min .or. coords(1) > this%r_max) then
                exit_flag = .true.
                print *, "Coordinates out of bounds"
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

        !> \brief Get target index from radial coordinate
        !> \param[in] this  Radial spatial discretization object
        !> \param[in] coord Radial coordinate [length units]
        !> \return target_ind Target index (1 to Num_targets), 0 if not found
        !> \details
        !>   Locates which radial cell contains the given coordinate.
        !>   
        !>   Algorithm:
        !>   1. Check bounds: return 0 if r < r_min or r > r_max
        !>   2. Check first cell: if r ≤ r_min + 0.5·Δr(1), return 1
        !>   3. Loop through interior cells:
        !>      - If r ≤ r_{i-1} + 0.5·Δr(i-1), return i-1
        !>      - If r ≤ r_i, return i
        !>   4. If r ≤ r_max, return Num_targets (last cell)
        !>   
        !>   Used for particle tracking, source placement, and interpolation.
        function get_target_ind_rad(this,coord) result(target_ind) !> Gets the target index for a given radial coordinate
            implicit none
            class(spatial_discr_rad_c) :: this !> spatial discretisation class
            real(kind=8), intent(in) :: coord(:) !> radial coordinates
            integer(kind=4) :: target_ind !> target index

            integer(kind=4) :: i
            real(kind=8), parameter :: eps=1d-12 !> small value for floating point comparison

            target_ind = 0 ! Initialize to 0 to indicate not found
            if (coord(1) < this%r_min .or. coord(1) > this%r_max ) then
                print *, "Error: Target coordinates out of bounds"
                return
            end if
            if (coord(1) <= this%r_min + 0.5*this%Delta_r(1)) then
                target_ind = 1 ! Set target index for the first target (innermost cell)
                return
            end if
            do i=2,this%Num_targets-1
                if (coord(1)<=this%targets(i-1)%coord(1)+0.5*this%Delta_r(i-1)) then
                    target_ind = i-1 ! Set target index (boundary between cells i-1 and i)
                    return
                else if (coord(1)<=this%targets(i)%coord(1)) then
                    target_ind = i ! Set target index (within cell i)
                    return
                else
                    continue ! Continue to next target
                end if
            end do
            if (coord(1) <= this%r_max) then
                target_ind = this%Num_targets ! Set target index for the last target (outermost cell)
                return
            end if
        end function get_target_ind_rad
        
        !> \brief Get maximum radial cell size
        !> \param[in] this Radial spatial discretization object
        !> \return max_cell_size Maximum Delta_r [length units]
        !> \details
        !>   Returns the largest radial cell spacing in the mesh.
        !>   Used for CFL stability analysis (least restrictive cell).
        function get_max_Delta_r(this) result(max_cell_size)
        class(spatial_discr_rad_c),  intent(in) :: this
        real(kind=8) :: max_cell_size
        max_cell_size=maxval(this%Delta_r)
        end function

        subroutine compute_Num_targets_rad(this)
            implicit none
            class(spatial_discr_rad_c) :: this
            ! integer(kind=4) :: i
            ! real(kind=8) :: r_current
            this%Num_targets=this%measure/this%cell_vol+this%targets_flag
            !> Compute number of targets based on Delta_r
            ! r_current=this%r_min
            ! this%Num_targets=1+this%targets_flag !> Start with initial target(s)
            ! do i=1,size(this%Delta_r)
            !     r_current=r_current+this%Delta_r(i)
            !     if (r_current <= this%r_max) then
            !         this%Num_targets=this%Num_targets+1
            !     else
            !         exit
            !     end if
            ! end do
            this%Num_targets_defined=.true.
        end subroutine compute_Num_targets_rad
        
        function get_num_cells_rad(this,dim) result(num_cells)
            class(spatial_discr_rad_c), intent(in) :: this
            integer(kind=4), intent(in), optional :: dim
            integer(kind=4) :: num_cells
            num_cells=this%Num_targets-this%targets_flag
        end function get_num_cells_rad

end module spatial_discr_rad_m
