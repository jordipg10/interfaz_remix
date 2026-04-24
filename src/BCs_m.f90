!> @file BCs_m.f90
!> @brief Boundary conditions module for reactive transport modeling
!> @details This module defines the boundary conditions type and associated procedures
!>          for handling various boundary condition types in partial differential equations

!> @brief Boundary conditions module
!> @details Contains the BCs_1D_c type and related procedures for managing boundary conditions
!>          in transport and flow simulations
module BCs_m
    use spatial_discr_rad_m, only: spatial_discr_rad_c !< Import radial spatial discretization class from radial spatial discretization module
    use spatial_discr_m, only: spatial_discr_c !< Import spatial discretization class from spatial discretization module
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c !< Import 2D homogeneous Eulerian mesh class from 2D spatial discretization module
    use time_fct_m, only: time_fct_real_c !< Import time-dependent function class from time function module
    implicit none !< Enforce explicit variable declaration
    save !< Preserve module variables between calls
    private !< Make all module entities private by default
    
    type, public, abstract :: BCs_c
    !> @brief Boundary conditions derived type
    !> @details Stores all boundary condition information including types, values, and flags
    !>          for both inflow and outflow boundaries
        integer(kind=4), allocatable :: labels(:)         !< Boundary condition type labels (1:inflow, 2:outflow in 1D)
                                                !< - 1: Dirichlet (prescribed value)
                                                !< - 2: Neumann homogeneous (zero dispersive flux)
                                                !< - 3: Robin (prescribed mass flux)
    end type

    type, public, extends(BCs_c) :: BCs_1D_c
    !> @brief One-dimensional boundary conditions derived type
        logical :: evap                         !< Evaporation flag: .true. if evaporation is active
        real(kind=8) :: conc_inf                !< Concentration at inflow boundary [M/L³]
        real(kind=8) :: conc_out                !< Concentration at outflow boundary [M/L³]
        real(kind=8) :: flux_inf                !< Flux at inflow boundary [M/L²/T]
        real(kind=8) :: flux_out                !< Flux at outflow boundary [M/L²/T]
        real(kind=8) :: head_inf                !< Hydraulic head at inflow boundary [L]
        real(kind=8) :: head_inf_D              !< Dimensionless hydraulic head at inflow boundary [-]
        real(kind=8) :: head_out                !< Hydraulic head at outflow boundary [L]
        real(kind=8) :: head_out_D              !< Dimensionless hydraulic head at outflow boundary [-]
        real(kind=8) :: caudal_inf              !< Flow rate (discharge) at inflow boundary [L³/T]
        type(time_fct_real_c) :: flow_inf       !< Time-dependent flow rate function at inflow
    contains
        procedure :: set_labels_1D              !< Set boundary condition type labels
        procedure :: allocate_labels=>allocate_labels_1D         !< Allocate labels array
        procedure :: set_evap                   !< Set evaporation flag
        procedure :: set_flow_inf               !< Set time-dependent inflow function
        procedure :: read_BCs=>read_BCs_1D       !< Read boundary condition types from file
        procedure :: read_Dirichlet_BCs_conc    !< Read Dirichlet concentration boundary conditions
        procedure :: read_Dirichlet_BCs_head    !< Read Dirichlet head boundary conditions
        procedure :: read_Robin_BC_inflow       !< Read Robin boundary condition at inflow
        procedure :: read_caudal_inf            !< Read constant inflow discharge from file
        procedure :: set_conc_boundary          !< Set concentration boundary values
        procedure :: set_cst_flux_boundary      !< Set constant flux boundary values
        procedure :: compute_dimless_BCs        !< Compute dimensionless boundary conditions
        procedure :: compute_flux_inf=>compute_flux_inf_1D         !< Compute inflow flux
    end type
    
    type, public, extends(BCs_1D_c) :: BCs_2D_c
    !> @brief Two-dimensional boundary conditions derived type
        !logical :: evap                         !< Evaporation flag: .true. if evaporation is active
        !real(kind=8) :: conc_up                !< Concentration at upper boundary [M/L³]
        !real(kind=8) :: conc_down                !< Concentration at lower boundary [M/L³]
        real(kind=8) :: flux_inf_y                !< Flux in y direction at inflow boundary [L³/L²/T]
        real(kind=8) :: flux_out_y                !< Flux in y direction at outflow boundary [L³/L²/T]
    contains
        procedure :: set_labels_2D              !< Set boundary condition type labels
        procedure :: allocate_labels=>allocate_labels_2D         !< Allocate labels array
        !procedure :: set_evap                   !< Set evaporation flag
        procedure :: read_BCs=>read_BCs_2D       !< Read boundary condition types from file
        !procedure :: read_Dirichlet_BCs_conc    !< Read Dirichlet concentration boundary conditions
        !procedure :: set_conc_boundary          !< Set concentration boundary values
        !procedure :: set_cst_flux_boundary      !< Set constant flux boundary values
        procedure :: compute_flux_inf=>compute_flux_inf_2D         !< Compute inflow flux
    end type
    contains
        !> @brief Set boundary condition type labels
        !> @details Assigns the boundary condition type for inflow and outflow boundaries
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] BCs Integer array of size 2 containing BC type labels (1:inflow, 2:outflow)
        !> @note Valid BC types are 1 (Dirichlet), 2 (Neumann), and 3 (Robin)
        !> @warning Will stop execution with error if BC type > 3
        subroutine set_labels_1D(this,labels)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            integer(kind=4), intent(in) :: labels(:) !< Array with BC type labels for inflow and outflow
            if (size(labels) /= 2) then
                error stop "BC labels array must have size 2" !< Check if labels array size is 2
            else if (labels(1)>3 .or. labels(2)>3) then
                error stop "BCs not implemented yet" !< Check if BC types are valid (≤3)
            end if
            this%labels=labels !< Assign BC labels to object
        end subroutine set_labels_1D

        subroutine set_labels_2D(this,labels)
            implicit none !< Enforce explicit variable declaration
            class(BCs_2D_c) :: this !< BCs_2D_c object instance
            integer(kind=4), intent(in) :: labels(:) !< Array with BC type labels for boundaries
            if (size(labels) /= 4) then
                error stop "BC labels array must have size 4" !< Check if labels array size is 4
            else if (any(labels>3) .or. any(labels<1)) then
                error stop "BCs not implemented yet" !< Check if BC types are valid (≤3)
            end if
            this%labels=labels !< Assign BC labels to object
        end subroutine set_labels_2D

        subroutine allocate_labels_1D(this)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            allocate(this%labels(2)) !< Allocate labels array of size 2
        end subroutine

        subroutine allocate_labels_2D(this)
            implicit none !< Enforce explicit variable declaration
            class(BCs_2D_c) :: this !< BCs_2D_c object instance
            allocate(this%labels(4)) !< Allocate labels array of size 4
        end subroutine
        
        !> @brief Set evaporation flag
        !> @details Enables or disables evaporation in the boundary condition
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] evap Logical flag: .true. to enable evaporation, .false. to disable
        subroutine set_evap(this,evap)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            logical, intent(in) :: evap !< Evaporation flag to be set
            this%evap=evap !< Assign evaporation flag to object
        end subroutine 
        
        !> @brief Set concentration boundary conditions
        !> @details Assigns concentration values at inflow and outflow boundaries
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] conc_inf Concentration at inflow boundary [M/L³]
        !> @param[in] conc_out Concentration at outflow boundary [M/L³]
        subroutine set_conc_boundary(this,conc_inf,conc_out)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            real(kind=8), intent(in) :: conc_inf,conc_out !< Concentrations at inflow and outflow
            this%conc_inf=conc_inf !< Assign inflow concentration
            this%conc_out=conc_out !< Assign outflow concentration
        end subroutine
        
        !> @brief Set constant flux boundary conditions
        !> @details Assigns the same constant flux value to both inflow and outflow boundaries
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] flux Flux value [M/L²/T] to be applied at both boundaries
        subroutine set_cst_flux_boundary(this,flux)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            real(kind=8), intent(in) :: flux !< Flux value for boundaries
            this%flux_inf=flux !< Assign flux to inflow boundary
            this%flux_out=this%flux_inf !< Assign same flux to outflow boundary
        end subroutine
        
        !> @brief Read boundary condition types from file
        !> @details Reads BC type labels and evaporation flag from specified file
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] filename Name of file containing BC information
        !> @note File format: line 1 = BC labels (2 integers), line 2 = evaporation flag (logical)
        subroutine read_BCs_1D(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file

            integer :: iostat !< I/O status variable
            character(len=100) :: line !< Line buffer for reading

            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            do
                read(1,*) line
                !print *, line
                !if (iostat /= 0) exit
                if (trim(line) == 'BOUNDARY CONDITIONS') then
                    call this%allocate_labels() !< Allocate labels array
                    !open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
                    read(1,*) this%labels !< Read BC type labels from first line
                    read(1,*) this%evap !< Read evaporation flag from second line
                else if (trim(line) == 'end') then
                    exit                                        !< End of boundary conditions block
                else
                    continue
                end if
            end do
            close(1) !< Close input file
        end subroutine
        
        !> @brief Read Dirichlet concentration boundary conditions
        !> @details Reads concentration values at inflow and outflow boundaries from file
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] filename Name of file containing concentration BC data
        !> @note File format: line 1 = inflow concentration, line 2 = outflow concentration
        subroutine read_Dirichlet_BCs_conc(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file
            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            read(1,*) this%conc_inf !< Read inflow concentration from first line
            read(1,*) this%conc_out !< Read outflow concentration from second line
            close(1) !< Close input file
        end subroutine
        
        !> @brief Read Robin boundary condition at inflow
        !> @details Reads discharge (caudal) value for Robin BC at inflow boundary
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] filename Name of file containing Robin BC data
        !> @note File format: line 1 = discharge at inflow [L³/T]
        !> @note Flux and concentration lines are commented out in current implementation
        subroutine read_Robin_BC_inflow(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file
            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            !read(1,*) this%flux_inf !< [Commented] Read flux at inflow
            read(1,*) this%caudal_inf !< Read discharge (caudal) at inflow from first line
            !read(1,*) this%conc_inf !< [Commented] Read concentration at inflow (not used in WMA)
            close(1) !< Close input file
        end subroutine
        
        !> @brief Read constant inflow discharge from file
        !> @details Reads a constant discharge value at the inflow boundary
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] filename Name of file containing discharge data
        !> @note File format: line 1 = discharge value [L³/T]
        !> @note The conversion to flux is commented out in current implementation
        subroutine read_caudal_inf(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file

            character(len=100) :: label !< Line buffer for reading
            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            do
                read(1,*) label
                if (trim(label) == 'DISCHARGE INFLOW') then
                    read(1,*) this%caudal_inf !< Read discharge from first line
                else if (trim(label) == 'end') then
                    exit                                        !< End of discharge inflow block
                else
                    continue
                end if
            end do
            close(1) !< Close input file
            !this%flux_inf=this%caudal_inf !< [Commented] Convert discharge to flux
        end subroutine
        
        !> @brief Set time-dependent inflow function
        !> @details Assigns a time-varying flow function to the inflow boundary
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] flow_inf Time-dependent flow function object
        subroutine set_flow_inf(this,flow_inf)
        implicit none !< Enforce explicit variable declaration
        class(BCs_1D_c) :: this !< BCs_1D_c object instance
        type(time_fct_real_c), intent(in) :: flow_inf !< Time-dependent flow function
        this%flow_inf=flow_inf !< Assign flow function to object
        end subroutine
        
        !> @brief Read Dirichlet head boundary conditions
        !> @details Reads hydraulic head values at inflow and outflow boundaries from file
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] filename Name of file containing head BC data
        !> @note File format: line 1 = head at inflow [L], line 2 = head at outflow [L]
        subroutine read_Dirichlet_BCs_head(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file
            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            read(1,*) this%head_inf !< Read hydraulic head at inflow from first line
            read(1,*) this%head_out !< Read hydraulic head at outflow from second line
            close(1) !< Close input file
        end subroutine
        
        !> @brief Compute dimensionless boundary conditions
        !> @details Normalizes hydraulic head values by characteristic head to obtain dimensionless heads
        !> @param[in,out] this BCs_1D_c object instance
        !> @param[in] h_c Characteristic head [L] used for normalization
        !> @note Dimensionless head = actual head / characteristic head
        subroutine compute_dimless_BCs(this,h_c)
            implicit none !< Enforce explicit variable declaration
            class(BCs_1D_c) :: this !< BCs_1D_c object instance
            real(kind=8), intent(in) :: h_c !< Characteristic head for normalization [L]
            this%head_inf_D=this%head_inf/h_c !< Compute dimensionless head at inflow
            this%head_out_D=this%head_out/h_c !< Compute dimensionless head at outflow
        end subroutine

        !> \brief Compute inflow flux from flow rate
        !> \details
        !>   Converts volumetric flow rate (caudal_inf) to flux based on geometry:
        !>   
        !>   For 2D radial coordinates:
        !>   \f[
        !>   \text{flux} = \frac{Q}{2\pi r_{min} b}
        !>   \f]
        !>   Where:
        !>   - Q: volumetric flow rate
        !>   - r_min: inner radius
        !>   - b: aquifer thickness
        !>   
        !>   For Cartesian 1D:
        !>   \f[
        !>   \text{flux} = Q
        !>   \f]
        !>   (assumes unit cross-sectional area or Q already as specific discharge)
        !>   
        !>   Polymorphic dispatch based on mesh type.
        !>
        !> \param[inout] this PDE object
        subroutine compute_flux_inf_1D(this,mesh)
        implicit none
        class(BCs_1D_c) :: this                                               !< PDE object
        class(spatial_discr_c), intent(in) :: mesh                              !< Spatial discretization object
        real(kind=8), parameter :: pi=4d0*atan(1d0)                           !< π constant
        select type (mesh)                                !< Polymorphic dispatch
        class is (spatial_discr_rad_c)                                        !< Radial coordinates
            if (mesh%dim==2) then
                !> 2D radial: flux = Q / (2πr*b)
                this%flux_inf=this%caudal_inf/(2d0*pi*mesh%r_min*mesh%targets(1)%thickness) !> Specific discharge [L/T]
            end if
        class default                                                         !< Cartesian or other
            this%flux_inf=this%caudal_inf                             !< flux = Q (assumes unit area)
        end select
        end subroutine compute_flux_inf_1D
        
        subroutine compute_flux_inf_2D(this,mesh)
        implicit none
        class(BCs_2D_c) :: this                                               !< PDE object
        class(spatial_discr_c), intent(in) :: mesh                              !< Spatial discretization object
        select type (mesh)
        type is (mesh_2D_Euler_homog_c)
            !> 2D Cartesian: flux = Q (assumes unit area or Q already as specific discharge)
            this%flux_inf=this%caudal_inf !< Assign flux in y direction (implementation depends on geometry)
            this%flux_inf_y=0d0 !< [Placeholder] Assign flux in y direction (implementation depends on geometry)
        end select  
        end subroutine compute_flux_inf_2D

        !> @brief Read boundary condition types from file
        !> @details Reads BC type labels and evaporation flag from specified file
        !> @param[in,out] this BCs_2D_c object instance
        !> @param[in] filename Name of file containing BC information
        !> @note File format: line 1 = BC labels (4 integers), line 2 = evaporation flag (logical)
        subroutine read_BCs_2D(this,filename)
            implicit none !< Enforce explicit variable declaration
            class(BCs_2D_c) :: this !< BCs_2D_c object instance
            character(len=*), intent(in) :: filename !< Name of input file

            integer :: iostat !< I/O status variable
            character(len=100) :: line !< Line buffer for reading

            open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
            do
                read(1,*) line
                !print *, line
                !if (iostat /= 0) exit
                if (trim(line) == 'BOUNDARY CONDITIONS') then
                    call this%allocate_labels() !< Allocate labels array
                    !open(unit=1,file=filename,status='old',action='read') !< Open file for reading (must exist)
                    read(1,*) this%labels !< Read BC type labels from first line
                    read(1,*) this%evap !< Read evaporation flag from second line
                else if (trim(line) == 'end') then
                    exit                                        !< End of boundary conditions block
                else
                    continue
                end if
            end do
            close(1) !< Close input file
        end subroutine
end module