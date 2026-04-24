!> \file diff_props_heterog_m.f90
!> \brief Heterogeneous diffusion properties module
!> \details Defines the heterogeneous diffusion properties class for 1D diffusion and transport problems
!> with spatially variable properties. This module extends the base properties class to handle:
!> - Spatially variable porosity
!> - Spatially variable dispersion coefficients (at cell centers and interfaces)
!> - Source/sink terms
!>
!> Heterogeneous properties allow for realistic representation of natural subsurface systems where
!> porosity and diffusion/dispersion vary in space due to geological heterogeneity, layering, or
!> variable sediment properties.

module diff_props_heterog_m
    use properties_m, only: props_c !> Import base properties class
    use spatial_discr_m, only: spatial_discr_c !> Import spatial discretization class
    implicit none !> Enforce explicit variable declarations
    save !> Preserve module variables between procedure calls
    private !> Private module scope: internal details hidden from outside modules
    !> \class diff_props_heterog_1D_c
    !> \brief Heterogeneous 1D diffusion properties subclass
    !> \details Extends the base properties class (props_c) to store spatially variable diffusion properties.
    !> This class is used when porosity and/or dispersion coefficients vary in space, requiring arrays
    !> to store values at each computational cell. Key features:
    !> - Porosity array: φ(x) at cell centers
    !> - Dispersion at cell interfaces: D(x+Δx/2) for flux calculations
    !> - Dispersion at cell centers: D(x) for storage and analysis
    !>
    !> Inherits from props_c and adds spatial variability capabilities.
    
    type, public, extends(props_c) :: diff_props_heterog_1D_c !> Heterogeneous 1D diffusion properties subclass extending base properties class
        !> \var porosity
        !> \brief Porosity at cell centers
        !> \details Array storing porosity values φ at each computational cell center.
        !> Porosity represents the fraction of void space available for fluid and solute storage.
        !> Range: 0 < φ ≤ 1 (dimensionless, but often expressed as fraction or percentage).
        !> Units: [-]
        real(kind=8), allocatable :: porosity(:)   !> Porosity array at cell centers: φ(x_i) for i=1..N_cells [-]
        
        !> \var diffusion_int
        !> \brief Dispersion coefficient at cell interfaces
        !> \details Array storing dispersion/diffusion coefficients at cell interfaces (boundaries between cells).
        !> Used for computing diffusive/dispersive fluxes between adjacent cells. For a grid with N cells,
        !> there are N+1 interfaces (including domain boundaries). Units: [L²/T]
        real(kind=8), allocatable :: diff_int(:) !> Dispersion coefficient at cell interfaces: D(x_i+1/2) for flux calculations [L²/T]
        
        !> \var diffusion_cent
        !> \brief Dispersion coefficient at cell centers
        !> \details Array storing dispersion/diffusion coefficients at cell centers.
        !> Used for storage, analysis, and conversion to interface values. Units: [L²/T]
        real(kind=8), allocatable :: diff_cent(:) !> Dispersion coefficient at cell centers: D(x_i) [L²/T]
        
        ! real(kind=8) :: long_disp !> Commented out: longitudinal dispersivity (alpha_L) - would be used for velocity-dependent dispersion [L]
    contains
        !> \brief Set heterogeneous diffusion properties
        !> \details Assigns values to porosity and dispersion arrays
        procedure :: set_props_diff_heterog
        
        !> \brief Read heterogeneous diffusion properties from file
        !> \details Reads porosity, dispersion, and source term from input file
        procedure :: read_props=>read_props_diff_heterog_1D
        procedure :: allocate_diff_props=>allocate_diff_props_1D
        !> \brief Check if properties are homogeneous
        !> \details Determines if spatially variable properties are actually uniform
        !procedure :: are_props_homog=>are_diff_props_homog
    end type

    type, public, extends(diff_props_heterog_1D_c) :: diff_props_heterog_2D_c !> Heterogeneous radial diffusion properties subclass
        !real(kind=8), allocatable :: porosity(:,:)   !> Porosity array at cell centers: φ(r_i) for i=1..N_cells [-]
        real(kind=8), allocatable :: diff_int_y(:) !> Dispersion coefficient at cell interfaces: D(r_i+1/2) for flux calculations [L²/T]
        real(kind=8), allocatable :: diff_cent_y(:) !> Dispersion coefficient at cell centers: D(r_i) [L²/T]
    contains
        procedure :: read_props=>read_props_diff_heterog_2D
    end type
    
    contains
    
        !> \brief Set heterogeneous diffusion properties
        !> \details Assigns values to the porosity and dispersion arrays for heterogeneous diffusion problems.
        !> This subroutine performs validation to ensure array dimensions are consistent:
        !> - Porosity and diffusion_cent must have the same size (one value per cell)
        !> - diffusion_int typically has N+1 values for N cells (includes interfaces at boundaries)
        !>
        !> \param[in,out] this Heterogeneous diffusion properties object to modify
        !> \param[in] porosity Porosity array at cell centers [-]
        !> \param[in] diffusion_int Dispersion coefficient array at cell interfaces [L²/T]
        !> \param[in] diffusion_cent Dispersion coefficient array at cell centers [L²/T]
        
        subroutine set_props_diff_heterog(this,porosity,diffusion_int,diffusion_cent)
            implicit none !> Enforce explicit variable declarations
            class(diff_props_heterog_1D_c) :: this !> Heterogeneous diffusion properties object (modified in place)
            real(kind=8), intent(in) :: porosity(:) !> Input porosity array at cell centers: φ(x_i) [-]
            real(kind=8), intent(in) :: diffusion_int(:) !> Input dispersion array at cell interfaces: D(x_i+1/2) [L²/T]
            real(kind=8), intent(in) :: diffusion_cent(:) !> Input dispersion array at cell centers: D(x_i) [L²/T]
            if (size(porosity)/=size(diffusion_cent)) error stop "Dimensions of porosity and dispersion at cell centres & 
                must be the same" !> Check array size consistency: porosity and diffusion_cent must have same length (one per cell)
            this%porosity=porosity !> Assign input porosity array to object's porosity [-]
            this%diff_int=diffusion_int !> Assign input interface dispersion array to object's diffusion_int [L²/T]
            this%diff_cent=diffusion_cent !> Assign input centered dispersion array to object's diffusion_cent [L²/T]
        end subroutine !> End of set_props_diff_heterog subroutine
        
        !> \brief Read heterogeneous diffusion properties from input file
        !> \details Reads diffusion properties from a formatted data file with the naming convention:
        !> root_diff_props.dat. The file contains:
        !> - Source/sink term: either uniform (flag=1, single value) or spatially variable (flag/=1, array)
        !> - Porosity: either uniform (flag=1, single value) or spatially variable (flag/=1, array)
        !> - Dispersion at cell centers: either uniform (flag=1, single value) or spatially variable (flag/=1, array)
        !>
        !> File format:
        !> Line 1: (blank)
        !> Line 2: flag r (for source term: flag=1 means uniform value r, flag/=1 means read array)
        !> Line 3: Array of source terms (if flag/=1) or skipped (if flag=1)
        !> Line 4: flag (for porosity)
        !> Line 5: flag phi (if flag=1: uniform porosity phi)
        !> Line 6: flag (for dispersion)
        !> Line 7: flag D (if flag=1: uniform dispersion D)
        !>
        !> \param[in,out] this Heterogeneous diffusion properties object to populate
        !> \param[in] root Root name for input file (file will be root_diff_props.dat)
        !> \param[in] spatial_discr Spatial discretization object containing grid information (number of targets, flags)
        
        subroutine read_props_diff_heterog_1D(this,root,spatial_discr)
            implicit none !> Enforce explicit variable declarations
            class(diff_props_heterog_1D_c) :: this !> Heterogeneous diffusion properties object to populate
            character(len=*), intent(in) :: root !> Root name for input file (filename = root_diff_props.dat)
            class(spatial_discr_c), intent(in), optional :: spatial_discr !> Spatial discretization object with grid information
            
            real(kind=8), parameter :: epsilon=1d-12 !> Small tolerance for numerical comparisons (not currently used) [-]
            real(kind=8) :: phi !> Temporary storage for uniform porosity value [-]
            real(kind=8) :: D !> Temporary storage for uniform dispersion coefficient [L²/T]
            real(kind=8) :: r !> Temporary storage for uniform source/sink term [M/(L³·T)]
            integer(kind=4) :: flag !> Flag indicating uniform (1) or variable (other) properties [-]
            
            open(unit=1,file=root//'_diff_props.dat',status='old',action='read') !> Open input file for reading: root_diff_props.dat
            read(1,"(/,F10.2)") flag !> Skip first line (/) then read flag for source term (format: floating point, 10 chars, 2 decimals) [-]
            if (flag == 1) then !> Check if source term is uniform (flag=1)
                backspace(1) !> Move file pointer back one line to re-read with different format
                read(1,*) flag, r !> Read flag and uniform source term value r [M/(L³·T)]
                allocate(this%source_term(spatial_discr%Num_targets-spatial_discr%targets_flag)) !> Allocate source term array with size = number of computational cells
                this%source_term=r !> Set all source term values to uniform value r [M/(L³·T)]
                this%source_term_order=0 !> Set source term order to 0 (zero-order/constant source) [-]
            else !> Source term is spatially variable
                read(1,*) this%source_term !> Read array of source term values (one per cell) [M/(L³·T)]
                if (size(this%source_term)/=spatial_discr%Num_targets-spatial_discr%targets_flag) error stop "Dimension error in & 
                    source term" !> Check array size matches number of computational cells
            end if !> End source term reading block
            read(1,*) flag !> Read flag for porosity: 1=uniform, other=variable [-]
            if (flag == 1) then !> Check if porosity is uniform (flag=1)
                backspace(1) !> Move file pointer back one line to re-read with different format
                read(1,*) flag, phi !> Read flag and uniform porosity value phi [-]
                allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag)) !> Allocate porosity array with size = number of computational cells
                this%porosity=phi !> Set all porosity values to uniform value phi [-]
            end if !> End porosity reading block (if flag/=1, porosity read elsewhere or remains unallocated)
            read(1,*) flag !> Read flag for dispersion: 1=uniform, other=variable [-]
            if (flag == 1) then !> Check if dispersion is uniform (flag=1)
                backspace(1) !> Move file pointer back one line to re-read with different format
                read(1,*) flag, D !> Read flag and uniform dispersion coefficient D [L²/T]
                allocate(this%diff_cent(spatial_discr%Num_targets-spatial_discr%targets_flag)) !> Allocate dispersion array at cell centers (workaround: chapuza)
                this%diff_cent=D !> Set all dispersion values to uniform value D (workaround: chapuza) [L²/T]
            end if !> End dispersion reading block
            close(1) !> Close input file
        end subroutine !> End of read_props_diff_heterog subroutine
        
      !> \brief Check if properties are homogeneous (spatially uniform)
      !> \details Determines whether the heterogeneous property arrays actually contain uniform values.
      !> This is useful for optimization: if properties are uniform despite being stored as arrays,
      !> simpler homogeneous solution methods can be used. The subroutine compares all array values
      !> to the first value within a tolerance (eps). If all values match, homog_flag is set to true.
      !>
      !> Properties checked:
      !> - Dispersion at cell centers: D(x_i) for all i
      !> - Porosity: φ(x_i) for all i
      !>
      !> \param[in,out] this Heterogeneous diffusion properties object (homog_flag is modified)
      
      subroutine are_diff_props_homog(this)
            implicit none !> Enforce explicit variable declarations
            class(diff_props_heterog_1D_c) :: this !> Heterogeneous diffusion properties object (homog_flag will be set)
            
            integer(kind=4) :: i !> Loop counter for iterating over cells [-]
            real(kind=8), parameter :: eps=1d-12 !> Tolerance for comparing property values (numerical epsilon) [-]
            
            this%homog_flag=.true. !> Initialize homogeneity flag to true (assume properties are uniform) [-]
            do i=2,size(this%porosity) !> Loop over all cells starting from second cell (compare to first cell)
                if (abs(this%diff_cent(1)-this%diff_cent(i))>eps .or. abs(this%porosity(1)-this%porosity(i))>eps) then !> Check if dispersion or porosity at cell i differs from cell 1 by more than tolerance
                    this%homog_flag=.false. !> Properties are heterogeneous: set flag to false [-]
                    exit !> Exit loop early (no need to check remaining cells)
                end if !> End comparison block
            end do !> End loop over cells
        end subroutine !> End of are_diff_props_homog subroutine

        subroutine read_props_diff_heterog_2D(this,root,spatial_discr)
            implicit none
            class(diff_props_heterog_2D_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
            
            real(kind=8), parameter :: epsilon=1d-12
            real(kind=8) :: phi !> Temporary storage for uniform porosity value [-]
            real(kind=8) :: D !> Temporary storage for uniform dispersion coefficient [L²/T]
            real(kind=8) :: r !> Temporary storage for uniform source/sink term [M/(L³·T)]
            integer(kind=4) :: flag !> Flag indicating uniform (1) or variable (other) properties [-]
            
            open(unit=1,file=root//'_diff_props.dat',status='old',action='read')
            read(1,"(/,F10.2)") flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, r
                allocate(this%source_term(spatial_discr%Num_targets))
                this%source_term=r
                this%source_term_order=0
            else
                read(1,*) this%source_term
                if (size(this%source_term)/=spatial_discr%Num_targets) error stop "Dimension error in source term"
            end if
            read(1,*) flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, phi
                allocate(this%porosity(spatial_discr%Num_targets))
                this%porosity=phi
            end if
            read(1,*) flag
            if (flag == 1) then
                backspace(1)
                read(1,*) flag, D
                allocate(this%diff_cent(spatial_discr%Num_targets))
                this%diff_cent=D
            end if
            close(1)
        end subroutine
        
        subroutine allocate_diff_props_1D(this,Num_cells)
            implicit none
            class(diff_props_heterog_1D_c) :: this
            integer(kind=4), intent(in) :: Num_cells
            
            call this%allocate_props(Num_cells) !> Call base class allocation for common properties (source term, etc.)
            allocate(this%porosity(Num_cells))
            allocate(this%diff_int(Num_cells+1))
            allocate(this%diff_cent(Num_cells))
    end subroutine
        
        
end module !> End of diff_props_heterog_m module