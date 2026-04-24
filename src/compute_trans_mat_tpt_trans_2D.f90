!> \file compute_trans_mat_tpt_trans_2D.f90
!! \brief Compute transition matrix T for transient 1D transport using Eulerian formulation
!! \details This subroutine constructs the spatial operator (transition matrix) for the transient
!! advection-dispersion equation in 1D using finite differences on an Eulerian mesh.
!!
!! The discretized transport equation has the form:
!! \f[
!!   \mathbf{F} \frac{d\mathbf{c}}{dt} = \mathbf{T} \mathbf{c} + \mathbf{g}
!! \f]
!! where:
!! - \f$\mathbf{F}\f$ is the mass/accumulation matrix (diagonal with porosity) [-]
!! - \f$\mathbf{T}\f$ is the transition matrix (tridiagonal, negative semi-definite) [1/T]
!! - \f$\mathbf{c}\f$ is the concentration vector [C]
!! - \f$\mathbf{g}\f$ is the source/boundary term vector [C/T]
!!
!! Properties of transition matrix \f$\mathbf{T}\f$:
!! - Tridiagonal structure (subdiagonal, diagonal, superdiagonal)
!! - Negative semi-definite (maximum principle preservation)
!! - Row sums equal zero when source term r = 0 (mass conservation)
!!
!! The matrix entries depend on:
!! - Spatial discretization scheme (centered, upwind, Petchamé-Carrera 2024)
!! - Mesh type (homogeneous, heterogeneous, radial)
!! - Target location (cell-centered vs. edge-centered)
!! - Dimensionless vs. dimensional formulation
!!
!! Finite difference schemes:
!! - scheme = 1: Traditional centered finite differences (2nd order accurate)
!! - scheme = 2: Petchamé & Carrera (2024) interface-based scheme (improved mass conservation)
!! - scheme = 3: Upwind scheme (1st order, stable for advection-dominated flow)
!!
!! \param[in,out] this Transport object (transport_1D_transient_c) containing mesh, properties, and transition matrix
!!
!! \use transport_transient_m Provides transport_1D_transient_c type for transient 1D transport
!! \use spatial_discr_1D_m Provides mesh_1D_Euler_homog_c and mesh_1D_Euler_heterog_c types for 1D Eulerian meshes
!! \use spatial_discr_rad_m Provides spatial_discr_rad_c type for radial (cylindrical/spherical) discretization
!! \use transport_properties_heterog_m Provides heterogeneous transport properties (dispersion, flux)

subroutine compute_trans_mat_tpt_trans_2D(this)
!> T: transition matrix (tridiagonal, negative semi-definite) [1/T]
!> rows sum = 0 if r=0 (mass conservation)
!> F*dc/dt=T*c+g (semi-discrete transport equation)
    use transport_transient_m, only: transport_2D_transient_c !> Transient 1D transport class
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c !> 2D Eulerian mesh types (homogeneous/heterogeneous)
    use spatial_discr_rad_m, only: spatial_discr_rad_c !> Radial spatial discretization class (cylindrical/spherical symmetry)
    use transport_properties_heterog_m, only: tpt_props_heterog_2D_c !> Heterogeneous transport properties (dispersion, flux, source terms)
    implicit none !> Enforce explicit variable declaration (no implicit typing)
    
    class(transport_2D_transient_c) :: this !> Transport object containing mesh, properties, and transition matrix to be computed
    
    real(kind=8) :: sign_flux !> Sign of Darcy flux (+1 for positive flow, -1 for negative flow, used in upwind scheme) [-]
    real(kind=8) :: r_i_12 !> Radius at interface i+1/2 for radial symmetry (not currently used) [L]
    real(kind=8) :: a,b !> Auxiliary variables for non-uniform mesh spacing: a = Δxᵢ₋₁/₂, b = Δxᵢ₊₁/₂ [L]
    integer(kind=4) :: i,j,n,Num_cells !> Loop index i over cells [1..n], total number of cells n [-]
    
    n=this%spatial_discr%Num_targets !> Get number of target points (computational nodes) from mesh discretization [-]
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> Get total number of spatial nodes (cells or interfaces) from mesh discretization [-]

    if (this%spatial_discr%adapt_ref.eq.1) then !> Check if adaptive mesh refinement flag is enabled (adapt_ref = 1 means mesh has changed) [-]
        deallocate(this%trans_mat%sub, this%trans_mat%diag, this%trans_mat%super, &
                   this%trans_mat%sub2, this%trans_mat%super2)
        call this%allocate_trans_mat() !> Reallocate transition matrix arrays with new mesh size (n-1, n, n-1 entries)
    end if
    
    select type (mesh=>this%spatial_discr) !> Polymorphic type selection to determine mesh type (homogeneous/radial) and compute matrix accordingly
    type is (mesh_2D_Euler_homog_c) !> Case: 1D homogeneous Eulerian mesh with uniform spacing Δx
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then !> Dimensionless formulation, cell-centered targets (targets_flag=0)
            if (mesh%scheme<3 .and. mesh%scheme>0) then !> Centered finite differences (scheme 1 or 2) for dimensionless case
                this%trans_mat%sub=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal: Tᵢ,ᵢ₋₁ = 1/Δx² + 1/(2Δx) (diffusion + backward advection) [1/L² + 1/L]
                this%trans_mat%super=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal: Tᵢ,ᵢ₊₁ = 1/Δx² - 1/(2Δx) (diffusion + forward advection) [1/L² + 1/L]
            !else if (mesh%scheme.eq.2) then !> proposed by Petchamé & Carrera (2024) - alternative scheme (commented out)
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> First superdiagonal entry for boundary cell
            !    do i=2,n-1 !> Loop over interior cells
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal for interior cells
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal for interior cells
            !    end do
            !    this%trans_mat%sub(n-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Last subdiagonal entry for boundary cell
            else if (mesh%scheme.eq.3) then !> Upwind finite difference scheme (scheme=3) for stability in advection-dominated transport
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then !> Check if flux has consistent sign (all positive or all negative)
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1)) !> Extract flux sign: +1 if flux > 0 (forward flow), -1 if flux < 0 (backward flow) [-]
                    this%trans_mat%sub=1d0/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*1d0/mesh%Delta_x !> Subdiagonal with upwind bias: T_sub = 1/Δx² + (sgn+1)/2·1/Δx (adds advection for backward flow) [1/L²]
                    this%trans_mat%super=1d0/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*1d0/mesh%Delta_x !> Superdiagonal with upwind bias: T_super = 1/Δx² + (sgn-1)/2·1/Δx (adds advection for forward flow) [1/L²]
                end if
            else !> Invalid scheme number for dimensionless case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported scheme
            end if
        else if (mesh%targets_flag .eq. 0) then !> Dimensional formulation, cell-centered targets (targets_flag=0)
            if (mesh%scheme.eq.1) then !> Traditional centered finite differences (scheme=1) with dimensional transport properties
                !> sub/super (x-direction): size Ny*(Nx-1), flat index = (i-1)*(Nx-1)+j
                ! DEBUG: print variables used in sub/super x-direction loops
                print *, "[DEBUG] mesh%Num_cells_x =", mesh%Num_cells_x
                print *, "[DEBUG] mesh%Num_cells_y =", mesh%Num_cells_y
                print *, "[DEBUG] mesh%Delta_x     =", mesh%Delta_x
                print *, "[DEBUG] mesh%Delta_y     =", mesh%Delta_y
                ! DEBUG: trans_mat allocation status and sizes
                print *, "[DEBUG] trans_mat%num_cols               =", this%trans_mat%num_cols
                print *, "[DEBUG] allocated(trans_mat%diag)        =", allocated(this%trans_mat%diag)
                if (allocated(this%trans_mat%diag)) &
                    print *, "[DEBUG] size(trans_mat%diag)             =", size(this%trans_mat%diag), &
                             " | expected =", mesh%Num_cells_x*mesh%Num_cells_y
                print *, "[DEBUG] allocated(trans_mat%sub)         =", allocated(this%trans_mat%sub)
                if (allocated(this%trans_mat%sub)) &
                    print *, "[DEBUG] size(trans_mat%sub)              =", size(this%trans_mat%sub), &
                             " | expected =", mesh%Num_cells_y*(mesh%Num_cells_x-1)
                print *, "[DEBUG] allocated(trans_mat%super)       =", allocated(this%trans_mat%super)
                if (allocated(this%trans_mat%super)) &
                    print *, "[DEBUG] size(trans_mat%super)            =", size(this%trans_mat%super), &
                             " | expected =", mesh%Num_cells_y*(mesh%Num_cells_x-1)
                print *, "[DEBUG] allocated(trans_mat%sub2)        =", allocated(this%trans_mat%sub2)
                if (allocated(this%trans_mat%sub2)) &
                    print *, "[DEBUG] size(trans_mat%sub2)             =", size(this%trans_mat%sub2), &
                             " | expected =", (mesh%Num_cells_y-1)*mesh%Num_cells_x
                print *, "[DEBUG] allocated(trans_mat%super2)      =", allocated(this%trans_mat%super2)
                if (allocated(this%trans_mat%super2)) &
                    print *, "[DEBUG] size(trans_mat%super2)           =", size(this%trans_mat%super2), &
                             " | expected =", (mesh%Num_cells_y-1)*mesh%Num_cells_x
                ! DEBUG: input arrays
                print *, "[DEBUG] size(disp_tensor_int) =", &
                    size(this%tpt_props_heterog%disp_tensor_int,1), &
                    size(this%tpt_props_heterog%disp_tensor_int,2), &
                    size(this%tpt_props_heterog%disp_tensor_int,3)
                print *, "[DEBUG] disp_tensor_int(:,:,1) ="
                print *, this%tpt_props_heterog%disp_tensor_int(:,:,1)
                print *, "[DEBUG] size(flux_cent)  =", size(this%tpt_props_heterog%flux_cent)
                print *, "[DEBUG] flux_cent        =", this%tpt_props_heterog%flux_cent
                do i=1,mesh%Num_cells_y !> Loop over rows i = 1..Num_cells_y (all y-rows)
                    do j=1,mesh%Num_cells_x-1 !> Loop over x-interfaces j = 1..Num_cells_x-1
                        !> Padding added to keep the original comment intact
                        !> sub flat index: (i-1)*(Num_cells_x-1)+j couples cell (i,j+1) to its left neighbour (i,j)
                        !>
                        !> Dispersion term: disp_tensor_int(j+1, i, 1) = D_xx at x-interface (j+1) in row i
                        !>   - disp_tensor_int is allocated (Num_cells_x+1, Num_cells_y+1, 3)
                        !>   - x-interface index j+1 lies between columns j and j+1 (0-based: x = j*Delta_x)
                        !>   - y-index i selects the row (y-interface index equals row index for cell-centred targets)
                        !>   - component 1 = D_xx (longitudinal dispersion in x-direction)
                        !>   Contribution: D_xx(j+1,i) / Delta_x^2  [1/T]
                        !>
                        !> Advection term: flux_cent((i-1)*Num_cells_x + j+1) = q_x at centre of cell (i, j+1)
                        !>   - flux_cent is a 1D array of size Num_cells_x*Num_cells_y in row-major order
                        !>   - flat index (i-1)*Num_cells_x + (j+1) maps 2D cell (row i, col j+1) to 1D storage
                        !>   - centred scheme uses the flux at the destination cell centre (i, j+1) for the left-neighbour coupling
                        !>   Contribution: +q_x(i, j+1) / (2*Delta_x)  [1/T]  (positive: subdiagonal adds upstream dispersion+advection)
                        this%trans_mat%sub((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j+1,i,1)/(mesh%Delta_x**2) + & !> Dispersion: D_xx at x-interface j+1 in row i, divided by Δx² [1/T]
                            this%tpt_props_heterog%flux_cent(&
                            (i-1)*mesh%Num_cells_x+j+1)/(2*mesh%Delta_x) !> Advection: q_x at cell centre (i,j+1), divided by 2Δx [1/T]
                    end do
                end do
                do i=1,mesh%Num_cells_y !> Loop over rows i = 1..Num_cells_y (all y-rows)
                    do j=1,mesh%Num_cells_x-1 !> Loop over x-interfaces j = 1..Num_cells_x-1
                        !> super flat index: (i-1)*(Num_cells_x-1)+j couples cell (i,j) to its right neighbour (i,j+1)
                        !>
                        !> Dispersion term: disp_tensor_int(j+1, i, 1) = D_xx at x-interface (j+1) in row i
                        !>   - Same shared interface as the sub entry: between columns j and j+1
                        !>   - x-interface index j+1, y-index i, component 1 = D_xx (x-direction dispersion)
                        !>   Contribution: D_xx(j+1,i) / Delta_x^2  [1/T]
                        !>
                        !> Advection term: flux_cent((i-1)*Num_cells_x + j) = q_x at centre of cell (i, j)
                        !>   - flux_cent is a 1D array of size Num_cells_x*Num_cells_y in row-major order
                        !>   - flat index (i-1)*Num_cells_x + j maps 2D cell (row i, col j) to 1D storage
                        !>   - centred scheme uses the flux at the source cell centre (i, j) for the right-neighbour coupling
                        !>     (contrast with sub, which uses the destination cell (i, j+1))
                        !>   - sign is negative: forward flux carries solute away from cell (i,j) into (i,j+1)
                        !>   Contribution: -q_x(i, j) / (2*Delta_x)  [1/T]
                        this%trans_mat%super((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j+1,i,1)/(mesh%Delta_x**2) - & !> Dispersion: D_xx at x-interface j+1 in row i, divided by Δx² [1/T]
                            this%tpt_props_heterog%flux_cent(&
                            (i-1)*mesh%Num_cells_x+j)/(2*mesh%Delta_x) !> Advection: q_x at cell centre (i,j), divided by 2Δx, subtracted [1/T]
                    end do
                end do
                !> sub2/super2 (y-direction): size (Ny-1)*Nx, flat index = (i-1)*Nx+j
                do i=1,mesh%Num_cells_y-1 !> Loop over interior cells to compute diagonal entries
                    do j=1,mesh%Num_cells_x !> Loop over x-direction cells
                        this%trans_mat%sub2((i-1)*mesh%Num_cells_x+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i+1,2)/(mesh%Delta_y**2) + & !> Subdiagonal: T_j,j₋₁ = Dj/Δx² + q_j/(2Δx) (dispersion + advection) [1/T]
                            this%tpt_props_heterog%flux_cent_y(&
                            (i)*mesh%Num_cells_x+j)/(2*mesh%Delta_y) !> Add advection term using centered flux values [1/T]
                    end do
                end do
                do i=1,mesh%Num_cells_y-1 !> Loop over interior cells to compute diagonal entries
                    do j=1,mesh%Num_cells_x !> Loop over x-direction cells
                        this%trans_mat%super2((i-1)*mesh%Num_cells_x+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i+1,2)/(mesh%Delta_y**2) - & !> Superdiagonal: T_j,j₊₁ = Dj/Δx² - q_j/(2Δx) (dispersion - advection) [1/T]
                            this%tpt_props_heterog%flux_cent_y(&
                            (i-1)*mesh%Num_cells_x+j)/(2*mesh%Delta_y) !> Subtract advection term using centered flux values [1/T]
                    end do
                end do
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024) interface-based scheme (scheme=2) for improved mass balance
                this%trans_mat%super(1)=this%tpt_props_heterog%disp_int(2)/(mesh%Delta_x**2) - & !> First superdiagonal: use interface dispersion D₁₊₁/₂ at x = x₁₊₁/₂ [1/T]
                    this%tpt_props_heterog%flux_int(2)/(2*mesh%Delta_x) !> Subtract interface flux q₁₊₁/₂ at x = x₁₊₁/₂ [1/T]
                do i=2,n-1 !> Loop over interior cells i=2..n-1 to compute matrix entries using interface values
                    this%trans_mat%sub(i-1)=this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) + & !> Subdiagonal: use dispersion at left interface Dᵢ₋₁/₂ [1/T]
                        this%tpt_props_heterog%flux_int(i)/(2*mesh%Delta_x) !> Add flux at left interface qᵢ₋₁/₂ [1/T]
                    this%trans_mat%super(i)=this%tpt_props_heterog%disp_int(i+1)/(mesh%Delta_x**2) - & !> Superdiagonal: use dispersion at right interface Dᵢ₊₁/₂ [1/T]
                        this%tpt_props_heterog%flux_int(i+1)/(2*mesh%Delta_x) !> Subtract flux at right interface qᵢ₊₁/₂ [1/T]
                end do
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%disp_int(n)/(mesh%Delta_x**2) + & !> Last subdiagonal: use interface dispersion Dₙ₋₁/₂ at x = xₙ₋₁/₂ [1/T]
                    this%tpt_props_heterog%flux_int(n)/(2*mesh%Delta_x) !> Add interface flux qₙ₋₁/₂ at x = xₙ₋₁/₂ [1/T]
            else if (mesh%scheme.eq.3) then !> Upwind scheme (scheme=3) for dimensional case with advection dominance
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1))
                    this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*&
                        this%tpt_props_heterog%flux_cent(2:n)/mesh%Delta_x
                    this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*&
                        this%tpt_props_heterog%flux_cent(1:n-1)/mesh%Delta_x
                end if
            else !> Invalid scheme number for dimensional cell-centered case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported scheme
            end if
        else if (this%dimless.eqv..false. .and. this%tpt_props_heterog%homog_flag.eqv..true.) & !> Edge-centered targets, dimensional formulation, homogeneous properties
            then !> (targets located at cell interfaces rather than centers)
            this%trans_mat%sub=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2) + & !> Subdiagonal: Tᵢ,ᵢ₋₁ = D_int,i/Δx² + q_int,i/(2Δx) using interface values [1/T]
                this%tpt_props_heterog%flux_int(2:n)/(2*mesh%Delta_x) !> Add advection term with interface flux values [1/T]
            this%trans_mat%super=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2) - & !> Superdiagonal: Tᵢ,ᵢ₊₁ = D_int,i/Δx² - q_int,i/(2Δx) using interface values [1/T]
                this%tpt_props_heterog%flux_int(2:n)/(2*mesh%Delta_x) !> Subtract advection term with interface flux values [1/T]
        end if
        this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term !> Initialize diagonal with negative source term: Tᵢᵢ = -fᵢ·rᵢ (flag multiplies source) [1/T]
        !> --- interior nodes: rows 2..Num_cells_y-1, cols 2..Num_cells_x-1 ---
        !> sub/super index for node (i,j): (i-1)*(Num_cells_x-1) + (j-1)
        !> sub2/super2 index for node (i,j): (i-1)*Nx + j
        do i=2,mesh%Num_cells_y-1
            do j=2,mesh%Num_cells_x-1
                this%trans_mat%diag((i-1)*mesh%Num_cells_x+j)= &
                    this%trans_mat%diag((i-1)*mesh%Num_cells_x+j) &
                    - this%trans_mat%sub( (i-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                    - this%trans_mat%super( (i-1)*(mesh%Num_cells_x-1)+j ) &
                    - this%trans_mat%sub2( (i-2)*mesh%Num_cells_x+j ) &
                    - this%trans_mat%super2( (i-1)*mesh%Num_cells_x+j )
            end do
        end do
        !> --- top boundary row (i=1): cols 2..Num_cells_x-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat%diag(j)= &
                this%trans_mat%diag(j) &
                - this%trans_mat%sub( j-1 ) &
                - this%trans_mat%super( j ) &
                - this%trans_mat%super2( j )
        end do
        !> --- bottom boundary row (i=Num_cells_y): cols 2..Num_cells_x-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j)= &
                this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j) &
                - this%trans_mat%sub( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                - this%trans_mat%super( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+j ) &
                - this%trans_mat%sub2( (mesh%Num_cells_y-2)*mesh%Num_cells_x+j )
        end do

    end select !> End polymorphic mesh type selection
    this%trans_mat_prev%sub=this%trans_mat%sub !> Store current subdiagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%diag=this%trans_mat%diag !> Store current diagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%super=this%trans_mat%super
    this%trans_mat_prev%sub2=this%trans_mat%sub2 !> Store current subdiagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%diag=this%trans_mat%diag !> Store current diagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%super=this%trans_mat%super
    this%trans_mat_prev%diag=this%trans_mat%diag !> Store current source term in previous transition matrix for time-stepping comparisons
    this%trans_mat_prev%super=this%trans_mat%super !> Store current source term in previous transition matrix for time-stepping comparisons
    this%trans_mat_prev%super2=this%trans_mat%super2 !> Store current source term in previous transition matrix for time-stepping comparisons
end subroutine !> End of compute_trans_mat_tpt_trans_Euler subroutine end subroutine !> End of compute_trans_mat_tpt_trans_Euler subroutine