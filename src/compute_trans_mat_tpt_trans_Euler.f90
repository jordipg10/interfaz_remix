!> \file compute_trans_mat_tpt_trans_Euler.f90
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

subroutine compute_trans_mat_tpt_trans_Euler(this)
!> T: transition matrix (tridiagonal, negative semi-definite) [1/T]
!> rows sum = 0 if r=0 (mass conservation)
!> F*dc/dt=T*c+g (semi-discrete transport equation)
    use transport_transient_m, only: transport_1D_transient_c !> Transient 1D transport class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c !> 1D Eulerian mesh types (homogeneous/heterogeneous)
    use spatial_discr_rad_m, only: spatial_discr_rad_c !> Radial spatial discretization class (cylindrical/spherical symmetry)
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c !> Heterogeneous transport properties (dispersion, flux, source terms)
    implicit none !> Enforce explicit variable declaration (no implicit typing)
    
    class(transport_1D_transient_c) :: this !> Transport object containing mesh, properties, and transition matrix to be computed
    
    real(kind=8) :: sign_flux !> Sign of Darcy flux (+1 for positive flow, -1 for negative flow, used in upwind scheme) [-]
    real(kind=8) :: r_i_12 !> Radius at interface i+1/2 for radial symmetry (not currently used) [L]
    real(kind=8) :: a,b !> Auxiliary variables for non-uniform mesh spacing: a = Δxᵢ₋₁/₂, b = Δxᵢ₊₁/₂ [L]
    integer(kind=4) :: i,n !> Loop index i over cells [1..n], total number of cells n [-]
    
    n=this%spatial_discr%Num_targets !> Get total number of spatial nodes (cells or interfaces) from mesh discretization [-]

    if (this%spatial_discr%adapt_ref.eq.1) then !> Check if adaptive mesh refinement flag is enabled (adapt_ref = 1 means mesh has changed) [-]
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super) !> Deallocate old transition matrix arrays (sub, diag, super) to prepare for new mesh
        call this%allocate_trans_mat() !> Reallocate transition matrix arrays with new mesh size (n-1, n, n-1 entries)
    end if
    
    select type (mesh=>this%spatial_discr) !> Polymorphic type selection to determine mesh type (homogeneous/radial) and compute matrix accordingly
    type is (mesh_1D_Euler_homog_c) !> Case: 1D homogeneous Eulerian mesh with uniform spacing Δx
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
                this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2) + & !> Subdiagonal: Tᵢ,ᵢ₋₁ = Dᵢ/Δx² + qᵢ/(2Δx) (dispersion + advection) [1/T]
                    this%tpt_props_heterog%flux_cent(2:n)/(2*mesh%Delta_x) !> Add advection term using centered flux values [1/T]
                this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) - & !> Superdiagonal: Tᵢ,ᵢ₊₁ = Dᵢ/Δx² - qᵢ/(2Δx) (dispersion - advection) [1/T]
                    this%tpt_props_heterog%flux_cent(1:n-1)/(2*mesh%Delta_x) !> Subtract advection term using centered flux values [1/T]
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
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then !> Check if flux has uniform sign across domain (no flow reversal)
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1)) !> Determine flow direction: +1 for forward (left-to-right), -1 for backward (right-to-left) [-]
                    this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*& !> Subdiagonal with upwind advection: D/Δx² + (sgn+1)/2·q/Δx [1/T]
                    this%tpt_props_heterog%flux_cent(2:n)/mesh%Delta_x !> Add upwind advection contribution (backward flow uses upwind from i-1)
                    this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*& !> Superdiagonal with upwind advection: D/Δx² + (sgn-1)/2·q/Δx [1/T]
                    this%tpt_props_heterog%flux_cent(1:n-1)/mesh%Delta_x !> Add upwind advection contribution (forward flow uses upwind from i+1)
                end if
            else !> Invalid scheme number for dimensional cell-centered case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported scheme
            end if
        else if (this%dimless.eqv..false. .and. this%tpt_props_heterog%homog_flag.eqv..true.) & !> Edge-centered targets, dimensional formulation, homogeneous properties
            then
            this%trans_mat%sub=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2) + & !> Subdiagonal: Tᵢ,ᵢ₋₁ = D_int,i/Δx² + q_int,i/(2Δx) using interface values [1/T]
                this%tpt_props_heterog%flux_int(2:n)/(2*mesh%Delta_x) !> Add advection term with interface flux values [1/T]
            !this%trans_mat_prev%sub=this%trans_mat%sub
            this%trans_mat%super=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2) - & !> Superdiagonal: Tᵢ,ᵢ₊₁ = D_int,i/Δx² - q_int,i/(2Δx) using interface values [1/T]
                this%tpt_props_heterog%flux_int(2:n)/(2*mesh%Delta_x) !> Subtract advection term with interface flux values [1/T]
            !this%trans_mat_prev%super=this%trans_mat%super
        end if
    type is (spatial_discr_rad_c) !> Case: Radial spatial discretization (cylindrical or spherical symmetry) with radial coordinate r
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then !> Dimensionless formulation, cell-centered targets (not implemented for radial)
            !if (mesh%scheme.eq.1) then !> Traditional centered finite differences (commented out - not implemented)
            !    this%trans_mat%sub=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal entry (placeholder)
            !    this%trans_mat%super=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal entry (placeholder)
            !else if (mesh%scheme.eq.2) then !> Petchamé & Carrera scheme (commented out - not implemented)
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> First superdiagonal
            !    do i=2,n-1 !> Loop over interior cells
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal
            !    end do
            !    this%trans_mat%sub(n-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Last subdiagonal
            !else if (mesh%scheme.eq.3) then !> Upwind scheme (commented out - not implemented)
            !    if (minval(this%tpt_props_heterog%flux)>=0d0 .or. maxval(this%tpt_props_heterog%flux)<0d0) then !> Check flux sign
            !        sign_flux=sign(1d0,this%tpt_props_heterog%flux(1)) !> Extract flux direction
            !        this%trans_mat%sub=1d0/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*1d0/mesh%Delta_x !> Subdiagonal with upwind
            !        this%trans_mat%super=1d0/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*1d0/mesh%Delta_x !> Superdiagonal with upwind
            !    end if
            !else !> Invalid scheme
            !    error stop "Scheme not implemented yet" !> Error termination
            !end if
        else if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..false.) then !> Dimensional formulation, cell-centered targets for radial geometry
            if (mesh%scheme.eq.1) then !> Traditional centered finite differences (scheme=1) for radial case (commented out - not implemented)
                !this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2) + this%tpt_props_heterog%flux(2:n)/& !> Subdiagonal (placeholder)
                !(2*mesh%Delta_x) !> Continuation of subdiagonal
                !this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2) - & !> Superdiagonal (placeholder)
                !this%tpt_props_heterog%flux(1:n-1)/(2*mesh%Delta_x) !> Continuation of superdiagonal
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024) interface-based scheme for non-uniform radial mesh
                this%trans_mat%super(1)=8d0*this%tpt_props_heterog%disp_int(2)/((mesh%Delta_r(1)+mesh%Delta_r(2))*& !> First superdiagonal using interface dispersion and non-uniform spacing [1/T]
                    (3d0*mesh%Delta_r(1)+mesh%Delta_r(2))) - & !> Denominator: product of averaged spacing and total spacing for first cell
                    this%tpt_props_heterog%flux_int(2)/(mesh%Delta_r(1)+mesh%Delta_r(2)) !> Subtract advection contribution using interface flux [1/T]
                do i=2,n-1 !> Loop over interior radial cells i=2..n-1 to compute transition matrix entries
                    a=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i-1)) !> Compute half-spacing to left interface: a = (Δrᵢ + Δrᵢ₋₁)/2 [L]
                    b=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1)) !> Compute half-spacing to right interface: b = (Δrᵢ + Δrᵢ₊₁)/2 [L]
                    this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%disp_int(i)/(a*(a+b)) + & !> Subdiagonal: Tᵢ,ᵢ₋₁ = 2Dᵢ₋₁/₂/(a·(a+b)) + qᵢ₋₁/₂/(2a) [1/T]
                        this%tpt_props_heterog%flux_int(i)/(2*a) !> Add advection contribution normalized by left half-spacing [1/T]
                    this%trans_mat%super(i)=2d0*this%tpt_props_heterog%disp_int(i+1)/(b*(a+b)) - & !> Superdiagonal: Tᵢ,ᵢ₊₁ = 2Dᵢ₊₁/₂/(b·(a+b)) - qᵢ₊₁/₂/(2b) [1/T]
                        this%tpt_props_heterog%flux_int(i+1)/(2*b) !> Subtract advection contribution normalized by right half-spacing [1/T]
                end do
                this%trans_mat%sub(n-1)=8d0*this%tpt_props_heterog%disp_int(n)/((mesh%Delta_r(n-1)+mesh%Delta_r(n))*& !> Last subdiagonal using interface dispersion and non-uniform spacing [1/T]
                    (3d0*mesh%Delta_r(n)+mesh%Delta_r(n-1))) + & !> Denominator: product of averaged spacing and total spacing for last cell
                    this%tpt_props_heterog%flux_int(n)/(mesh%Delta_r(n-1)+mesh%Delta_r(n)) !> Add advection contribution using interface flux [1/T]
            else if (mesh%scheme.eq.3) then !> Upwind scheme (scheme=3) for radial case (commented out - not fully implemented)
                !if (minval(this%tpt_props_heterog%flux)>=0d0 .or. maxval(this%tpt_props_heterog%flux)<0d0) then !> Check for uniform flux direction
                !    sign_flux=sign(1d0,this%tpt_props_heterog%flux(1)) !> Extract flux sign
                !    this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*& !> Subdiagonal with upwind
                !    this%tpt_props_heterog%flux(2:n)/mesh%Delta_x !> Upwind advection term
                !    this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*& !> Superdiagonal with upwind
                !    this%tpt_props_heterog%flux(1:n-1)/mesh%Delta_x !> Upwind advection term
                !end if
            else !> Invalid scheme number for radial case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported radial scheme
            end if
        else if (mesh%targets_flag.eq.1 .and. this%dimless.eqv..false.) then !> Dimensional formulation, edge-centered targets for radial mesh (targets at interfaces)
            this%trans_mat%super(1)=this%tpt_props_heterog%disp_cent(1)/(mesh%Delta_r(1)**2) - & 
                this%tpt_props_heterog%flux_int(1)/(2d0*mesh%Delta_r(1)) !> First superdiagonal: T₁,₂ = 2D₁/Δr₁² for edge-centered radial discretization [1/T]
                do i=2,n-1 !> Loop over interior interface nodes i=2..n-1 for edge-centered radial Lagrangian formulation
                    !a=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i-1)) !> Half-spacing computation (commented out - not used in this formulation)
                    !b=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1)) !> Half-spacing computation (commented out - not used in this formulation)
                    this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%disp_cent(i-1)/( & !> Subdiagonal: Tᵢ,ᵢ₋₁ = 2Dᵢ / [Δrᵢ₋₁·(Δrᵢ₋₁+Δrᵢ)] (dispersion only) [1/T]
                        mesh%Delta_r(i-1)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) + &
                        this%tpt_props_heterog%flux_int(i)*mesh%Delta_r(i-1)/(mesh%Delta_r(i)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) !> Denominator: product of left cell size and total spacing
                    this%trans_mat%super(i)=2d0*this%tpt_props_heterog%disp_cent(i)/( & !> Superdiagonal: Tᵢ,ᵢ₊₁ = 2Dᵢ / [Δrᵢ·(Δrᵢ₋₁+Δrᵢ)] (dispersion only) [1/T]
                        mesh%Delta_r(i)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) - &
                        this%tpt_props_heterog%flux_int(i)*mesh%Delta_r(i)/(mesh%Delta_r(i-1)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) !> Denominator: product of right cell size and total spacing
                end do
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%disp_cent(n-1)/(mesh%Delta_r(n-1)**2) + &
                    this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_r(n-1)) !> Last subdiagonal: Tₙ,ₙ₋₁ = 2Dₙ/Δrₙ₋₁² for edge-centered radial discretization [1/T]
        else !> Invalid configuration for radial mesh (e.g., edge-centered targets or other unimplemented options)
            error stop "Radial transition matrix not fully implemented yet" !> Terminate program indicating radial implementation is incomplete
        end if
    end select !> End polymorphic mesh type selection
    this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term !> Initialize diagonal with negative source term: Tᵢᵢ = -fᵢ·rᵢ (flag multiplies source) [1/T]
    this%trans_mat%diag(2:n-1)=this%trans_mat%diag(2:n-1)-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1) !> Complete diagonal for interior nodes: Tᵢᵢ = -(Tᵢ,ᵢ₋₁ + Tᵢ,ᵢ₊₁) - r ensuring row sum = -r [1/T]
    !> Copy previous transition matrix
    this%trans_mat_prev%sub=this%trans_mat%sub !> Store current subdiagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%diag=this%trans_mat%diag !> Store current diagonal in previous transition matrix for time-stepping comparisons this%trans_mat_prev%super=this%trans_mat%super
    this%trans_mat_prev%diag=this%trans_mat%diag !> Store current source term in previous transition matrix for time-stepping comparisons
    this%trans_mat_prev%super=this%trans_mat%super !> Store current source term in previous transition matrix for time-stepping comparisons
end subroutine !> End of compute_trans_mat_tpt_trans_Euler subroutine 