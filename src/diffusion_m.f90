!****************************************************************************************************************************************************
!> \file diffusion_m.f90
!> \brief Module for 1D steady-state diffusion equation with source/sink terms
!>
!> \details
!> This module defines the `diffusion_1D_c` class for solving the one-dimensional steady-state
!> diffusion equation with spatially variable properties and source/sink terms.
!>
!> **Governing Equation:**
!>
!> The steady-state diffusion equation with source/sink term is:
!> \f[
!>    0 = D \frac{d^2c}{dx^2} + r(c_r - c)
!> \f]
!>
!> Where:
!> - \f$ c(x) \f$ = Concentration field [mol/L]
!> - \f$ D(x) \f$ = Diffusion coefficient [m²/s], may be spatially variable
!> - \f$ r(x) \f$ = Source/sink rate coefficient [1/s], spatially variable
!> - \f$ c_r(x) \f$ = Recharge/source concentration [mol/L]
!>
!> **Physical Interpretation:**
!>
!> - **Diffusion term**: \f$ D \nabla^2 c \f$ represents molecular or mechanical dispersion
!> - **Source/sink term**: \f$ r(c_r - c) \f$ represents:
!>   - When r > 0 and c_r > c: Source (e.g., recharge from surface)
!>   - When r > 0 and c_r < c: Sink (e.g., evapotranspiration, uptake)
!>   - When r = 0: Pure diffusion
!>
!> **Discretization:**
!>
!> The equation is discretized using finite differences on a 1D mesh, resulting in:
!> \f[
!>    (T + R) \mathbf{c} = \mathbf{g}
!> \f]
!>
!> Where:
!> - \f$ T \f$ = Transmission matrix (from diffusion term, tridiagonal)
!> - \f$ R \f$ = Recharge matrix (diagonal, R_{ii} = r_i)
!> - \f$ \mathbf{g} \f$ = Source term vector (includes boundary conditions)
!>
!> **Boundary Conditions:**
!>
!> - Dirichlet: c(0) = c_inf, c(L) = c_out (prescribed concentrations)
!> - Can handle heterogeneous properties and spatially variable recharge
!>
!> \see PDE_1D_c, diff_props_heterog_1D_c, solve_diff_1D, compute_trans_mat_diff
!****************************************************************************************************************************************************
module diffusion_m
    use PDE_m, only: PDE_1D_c                           !< Import base PDE class
    use diff_props_heterog_m, only: diff_props_heterog_1D_c
                                                        !< Import heterogeneous diffusion properties
    implicit none
    save    !< Preserve module variables between calls
    private !< Private module scope: internal details hidden from outside modules
!****************************************************************************************************************************************************
!> \brief Derived type for 1D steady-state diffusion equation
!>
!> \details
!> This type extends the base PDE_1D_c class to handle steady-state diffusion with:
!> - Spatially variable diffusion coefficients
!> - Source/sink terms with recharge concentrations
!> - Heterogeneous properties
!> - Boundary condition handling
!>
!> The class manages concentration fields, external concentrations (boundary + recharge),
!> and flags indicating where recharge occurs.
    type, public, extends(PDE_1D_c) :: diffusion_1D_c !< 1D diffusion equation class (extends PDE_1D_c)
    
        !> \subsection concentration_fields Concentration Fields
        real(kind=8), allocatable :: conc(:)        !< Solution concentration field [mol/L]
                                                    !< Dimension: [Num_targets]
                                                    !< Primary unknown being solved for
                                                    
        real(kind=8), allocatable :: conc_ext(:)    !< External concentrations [mol/L]
                                                    !< Dimension: [Num_targets]
                                                    !< Includes both recharge and boundary concentrations
                                                    !< Used in source term: g = R·c_ext
                                                    
        real(kind=8), allocatable :: conc_r(:)      !< Recharge concentrations [mol/L]
                                                    !< Dimension: [Num_targets]
                                                    !< c_r in source term r(c_r - c)
                                                    !< Only relevant where r > 0
                                                    
        real(kind=8), allocatable :: conc_bd(:)     !< Boundary concentrations [mol/L]
                                                    !< Dimension: [Num_targets]
                                                    !< Prescribed concentrations at boundaries
                                                    !< c_bd > 0 for Dirichlet boundaries
                                                    
        !> \subsection flags Recharge Flags
        integer(kind=4), allocatable :: conc_r_flag(:) !< Recharge indicator flag [-]
                                                    !< Dimension: [Num_targets]
                                                    !< = 1 if r > 0 (recharge/sink present)
                                                    !< = 0 otherwise (pure diffusion)
                                                    
        !> \subsection properties Heterogeneous Diffusion Properties
        type(diff_props_heterog_1D_c) :: diff_props_heterog !< Spatially variable diffusion properties
                                                    !< Contains: D(x), r(x), porosity, etc.
    contains
        !> \subsection setters Setter Procedures
        procedure :: set_conc_ext           !< Set external concentration field
        procedure :: set_diff_props_heterog !< Set heterogeneous diffusion properties
        procedure :: set_conc_r_flag=>set_conc_r_flag_diff
                                                    !< Set recharge flags based on r(x)
        
        !> \subsection computations Computational Procedures
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_diff
                                                    !< Compute transmission matrix T from diffusion
        procedure :: compute_source_term_PDE=>compute_source_term_diff
                                                    !< Compute source term g = R·c_ext + BCs
        procedure :: compute_rech_mat_PDE=>compute_rech_mat_diff
                                                    !< Compute recharge matrix R = diag(r)
        procedure :: initialise_PDE=>initialise_diffusion_1D
                                                    !< Initialize diffusion problem from input file
        procedure :: write_PDE=>write_diffusion_1D
                                                    !< Write solution output
        procedure :: allocate_conc          !< Allocate concentration arrays
        procedure :: solve_PDE=>solve_diff_1D
                                                    !< Solve steady-state diffusion equation
    end type
    
!****************************************************************************************************************************************************
!> \subsection interfaces Interface Definitions
!> Abstract interfaces for procedures implemented in separate files
    interface
        
        !> \brief Initialize 1D diffusion problem from input file
        !> \param[in,out] this Diffusion object to initialize
        !> \param[in] root Root name for input files (e.g., "problem" reads "problem_diff.dat")
        subroutine initialise_diffusion_1D(this,path,root,mesh_type)
            import diffusion_1D_c               !< Import type from parent module
            class(diffusion_1D_c) :: this       !< Diffusion object
            character(len=*), intent(in) :: path !> Input file path
            character(len=*), intent(in) :: root !> Root name for input files
            integer(kind=4), intent(in) :: mesh_type
        end subroutine
        
        !> \brief Compute transmission matrix from diffusion coefficients
        !> \param[in,out] this Diffusion object with updated transmission matrix
        subroutine compute_trans_mat_diff(this)
            import diffusion_1D_c               !< Import type from parent module
            implicit none
            class(diffusion_1D_c) :: this       !< Diffusion object
        end subroutine
        
        !> \brief Write diffusion solution to output file
        !> \param[in] this Diffusion object with solution
        !> \param[in] root Root name for output files
        !> \param[in] Time_out Time points for output [s]
        !> \param[in] output Solution matrix [Num_targets × Num_times]
       subroutine write_diffusion_1D(this)
            import diffusion_1D_c               !< Import type from parent module
            class(diffusion_1D_c), intent(in) :: this !< Diffusion object
            ! character(len=*), intent(in) :: root !> root name for output files
            ! real(kind=8), intent(in) :: Time_out(:) !< Output time points [s]
            ! real(kind=8), intent(in) :: output(:,:) !< Solution matrix [space × time]
        end subroutine
      
        !> \brief Solve 1D steady-state diffusion equation
        !> \param[in,out] this Diffusion object with problem data and solution
        !> \param[in] Time_out Time points for output [s] (single value for steady-state)
        !> \param[out] output Solution concentration field [mol/L]
        subroutine solve_diff_1D(this)
        import diffusion_1D_c               !< Import type from parent module
        class(diffusion_1D_c) :: this       !< Diffusion object
        !real(kind=8), intent(in) :: Time_out(:)  !< Output time points
        !real(kind=8), intent(out) :: output(:,:) !< Solution matrix
        end subroutine
        
    end interface
    
!****************************************************************************************************************************************************
!> \subsection implementations Module Procedure Implementations
    contains
        
!****************************************************************************************************************************************************
!> \brief Set external concentration field
!>
!> \details
!> Sets the external concentration field used in source/sink terms. The external concentrations
!> represent recharge concentrations where r > 0 and can include boundary concentrations.
!>
!> The source term is computed as: g = R·c_ext where R = diag(r)
!>
!> \param[in,out] this Diffusion object
!> \param[in] conc_ext External concentration array [mol/L], dimension must match Num_targets
!>
!> \warning Stops with error if dimension mismatch detected
        subroutine set_conc_ext(this,conc_ext)
            class(diffusion_1D_c) :: this           !< Diffusion object
            real(kind=8), intent(in) :: conc_ext(:) !< External concentrations [mol/L]
            
            !> Validate array dimensions
            if (size(conc_ext)/=this%spatial_discr%Num_targets) error stop "Dimension error in external concentration"
                                                    !< Check conc_ext size matches mesh size
            
            !> Assign external concentrations
            this%conc_ext=conc_ext                  !< Copy external concentration field
        end subroutine 
        
!****************************************************************************************************************************************************
!> \brief Set heterogeneous diffusion properties
!>
!> \details
!> Assigns spatially variable diffusion properties to the diffusion object, including:
!> - Diffusion coefficients D(x)
!> - Source/sink rates r(x)
!> - Porosity
!> - Other transport properties
!>
!> \param[in,out] this Diffusion object
!> \param[in] diff_props_heterog Heterogeneous diffusion properties object
        subroutine set_diff_props_heterog(this,diff_props_heterog)
            implicit none
            class(diffusion_1D_c) :: this           !< Diffusion object
            class(diff_props_heterog_1D_c), intent(in) :: diff_props_heterog
                                                    !< Heterogeneous property object
            
            !> Assign heterogeneous properties
            this%diff_props_heterog=diff_props_heterog !< Copy all property fields
        end subroutine
        
!****************************************************************************************************************************************************
!> \brief Set recharge flags based on source/sink rates
!>
!> \details
!> Initializes the recharge flag array to indicate where source/sink terms are active.
!> This is used to identify regions with recharge (r > 0) vs. pure diffusion (r = 0).
!>
!> Flag values:
!> - conc_r_flag(i) = 1 if r(i) > 0 (recharge/sink present)
!> - conc_r_flag(i) = 0 if r(i) = 0 (pure diffusion)
!>
!> \param[in,out] this Diffusion object with allocated conc_r_flag array
        subroutine set_conc_r_flag_diff(this)
            implicit none
            class(diffusion_1D_c) :: this           !< Diffusion object
            
            integer(kind=4) :: i                    !< Loop counter for grid points
            
            !> Allocate recharge flag array
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets))
                                                    !< Dimension: [Num_targets]
            
            !> Initialize all flags to zero (no recharge)
            this%conc_r_flag=0                      !< Default: pure diffusion everywhere
            
            !> Loop over all grid points and set flags
            do i=1,this%spatial_discr%Num_targets   !< Loop over mesh points
                !> Check if source/sink rate is positive
                if (this%diff_props_heterog%source_term(i)>0) then
                                                    !< If r(i) > 0: recharge or sink present
                    this%conc_r_flag(i)=1           !< Set flag indicating active source/sink
                end if
            end do                                  !< End loop over grid points
        end subroutine 
        
!****************************************************************************************************************************************************
!> \brief Compute source term vector for diffusion equation
!>
!> \details
!> Constructs the right-hand side source term vector g for the discrete diffusion equation:
!> (T + R)·c = g
!>
!> The source term includes:
!> 1. **Interior points**: g_i = r_i · c_ext,i (recharge contribution)
!> 2. **Boundary points**: 
!>    - g_1 = bd_mat(1) · c_inf (inflow boundary condition)
!>    - g_n = bd_mat(2) · c_out (outflow boundary condition)
!>
!> Mathematical formulation:
!> \f[
!>    g = R \cdot c_{ext} + g_{BC}
!> \f]
!>
!> Where:
!> - R = diag(r) is the recharge matrix
!> - c_ext = external (recharge) concentrations
!> - g_BC = boundary condition contributions
!>
!> \param[in,out] this Diffusion object with updated source_term_PDE vector
!>
!> \note Legacy code for alternative implementations is commented out
        subroutine compute_source_term_diff(this)
            !> Mathematical formula: g = R·c_ext + BC terms
            !use PDE_m, only: PDE_1D_c
            !use transport_m, only: transpoRT_c, diffusion_1D_c
            !use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
            !< (Commented) Legacy module imports
            implicit none
            class(diffusion_1D_c) :: this           !< Diffusion object
    
            !> Legacy allocation code (commented - array should be pre-allocated)
            !allocate(this%source_term_PDE(this%spatial_discr%Num_targets))
            !this%source_term_PDE=0d0 !> $g=0$ chapuza
            !< (Commented) Would allocate and initialize source term
            
            !> Legacy polymorphic selection (commented - now handled directly)
            !select type (this)
            !class is (diffusion_1D_c)
            !    this%source_term_PDE=this%rech_mat%diag*this%conc_ext
            !class is (diffusion_1D_transient_c)
            !    this%source_term_PDE=this%rech_mat%diag*this%conc_ext
            !end select
            !< (Commented) Type-specific handling
            
            !> \subsection interior_source Interior source term from recharge
            !> Compute element-wise product: g = R·c_ext where R = diag(r)
            this%source_term_PDE=this%rech_mat%diag*this%conc_ext
                                                    !< g_i = r_i · c_ext,i for all i
                                                    !< Element-wise multiplication
            
            !> \subsection boundary_sources Boundary condition contributions
            !> Inflow boundary (left/inner boundary)
            this%source_term_PDE(1)=this%bd_mat(1)*this%BCs%conc_inf
                                                    !< g_1 = bd_mat(1) · c_inf
                                                    !< Dirichlet BC at x=0
            
            !> Outflow boundary (right/outer boundary)
            this%source_term_PDE(this%spatial_discr%Num_targets)=this%bd_mat(2)*this%BCs%conc_out
                                                    !< g_n = bd_mat(2) · c_out
                                                    !< Dirichlet BC at x=L
            
            !> Legacy polymorphic code for different PDE types (commented)
            ! select type (this)
            ! type is (transpoRT_c)
            !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
            ! type is (transport_1D_transient_c)
            !     this%source_term_PDE=this%conc_r_flag*this%tpt_props_heterog%source_term*this%conc_ext
            ! type is (diffusion_1D_c)
            !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
            ! type is (diffusion_1D_transient_c)
            !     this%source_term_PDE=this%diff_props_heterog%source_term*this%conc_ext
            ! end select
            !< (Commented) Alternative implementations for derived types

        end subroutine
        
!****************************************************************************************************************************************************
!> \brief Allocate concentration arrays
!>
!> \details
!> Allocates all concentration-related arrays for the diffusion problem with proper dimensions
!> based on the spatial discretization (number of grid points).
!>
!> Arrays allocated:
!> - conc: Solution concentration field
!> - conc_ext: External (recharge) concentrations
!> - conc_r: Recharge concentrations
!> - conc_bd: Boundary concentrations
!> - conc_r_flag: Recharge indicator flags
!>
!> All arrays have dimension [Num_targets]
!>
!> \param[in,out] this Diffusion object with allocated concentration arrays
        subroutine allocate_conc(this)
            implicit none
            class(diffusion_1D_c) :: this           !< Diffusion object
            
            !> Allocate all concentration arrays in a single statement
            allocate(this%conc(this%spatial_discr%Num_targets), this%conc_ext(this%spatial_discr%Num_targets), &
                this%conc_r(this%spatial_discr%Num_targets), this%conc_bd(this%spatial_discr%Num_targets), &
                this%conc_r_flag(this%spatial_discr%Num_targets))
                                                    !< All arrays: dimension [Num_targets]
                                                    !< conc: solution field
                                                    !< conc_ext: external concentrations
                                                    !< conc_r: recharge concentrations
                                                    !< conc_bd: boundary concentrations
                                                    !< conc_r_flag: recharge flags
        end subroutine
        
!****************************************************************************************************************************************************
!> \brief Compute recharge matrix from source/sink rates
!>
!> \details
!> Constructs the recharge matrix R as a diagonal matrix with source/sink rates on the diagonal:
!> \f[
!>    R = \text{diag}(r_1, r_2, \ldots, r_n)
!> \f]
!>
!> Where r_i is the source/sink rate at grid point i [1/s].
!>
!> This matrix appears in the discrete diffusion equation:
!> \f[
!>    (T + R) \cdot c = g
!> \f]
!>
!> The R matrix contributes the source/sink term: r(c_r - c)
!>
!> \param[in,out] this Diffusion object with updated rech_mat diagonal
        subroutine compute_rech_mat_diff(this)
        implicit none
        class(diffusion_1D_c) :: this               !< Diffusion object
        
        !> Assign source/sink rates to recharge matrix diagonal
        !> Mathematical formula: R = diag(r)
        this%rech_mat%diag=this%diff_props_heterog%source_term
                                                    !< R_{ii} = r_i for i=1..n
                                                    !< Source/sink rate at each grid point
        end subroutine
        
end module !< End diffusion_m 