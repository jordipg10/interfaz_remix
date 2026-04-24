!> \file compute_trans_mat_tpt_trans_Lagr.f90
!! \brief Compute transition matrix T for transient 1D transport using Lagrangian formulation
!! \details This subroutine constructs the spatial operator (transition matrix) for the transient
!! transport equation in 1D using finite differences with Lagrangian framework (moving with flow).
!!
!! The Lagrangian formulation follows the fluid particles, removing the advection term from the
!! spatial operator. The discretized transport equation has the form:
!! \f[
!!   \mathbf{F} \frac{d\mathbf{c}}{dt} = \mathbf{T} \mathbf{c} + \mathbf{g}
!! \f]
!! where:
!! - \f$\mathbf{F}\f$ is the storage matrix (diagonal with porosity) [-]
!! - \f$\mathbf{T}\f$ is the transition matrix (tridiagonal, negative semi-definite) [1/T]
!! - \f$\mathbf{c}\f$ is the concentration vector [C]
!! - \f$\mathbf{g}\f$ is the external terms vector [C/T]
!!
!! **Key difference from Eulerian formulation:**
!! The Lagrangian approach tracks fluid particles moving with the flow, so the transition matrix
!! contains only dispersion terms (no advection). The advection is handled implicitly by the
!! moving coordinate system.
!!
!! Properties of transition matrix \f$\mathbf{T}\f$:
!! - Tridiagonal structure (subdiagonal, diagonal, superdiagonal)
!! - Negative semi-definite (maximum principle preservation)
!! - Row sums equal zero when source term r = 0 (mass conservation)
!! - Contains only dispersion/diffusion terms (no explicit advection)
!!
!! Finite difference schemes:
!! - scheme = 1: Traditional centered finite differences (2nd order accurate, dispersion only)
!! - scheme = 2: Petchamé & Carrera (2024) interface-based scheme (improved mass conservation)
!! - scheme = 3: Upwind scheme (dispersion only, advection handled by Lagrangian framework)
!!
!! Supported mesh types:
!! - 1D homogeneous Eulerian mesh (uniform spacing)
!! - 1D heterogeneous Eulerian mesh (variable spacing)
!! - Radial mesh (cylindrical/spherical symmetry with non-uniform spacing)
!!
!! \param[in,out] this Transport object (transport_1D_transient_c) containing mesh, properties, and transition matrix
!!
!! \use transport_transient_m Provides transport_1D_transient_c type for transient 1D transport
!! \use spatial_discr_1D_m Provides mesh_1D_Euler_homog_c and mesh_1D_Euler_heterog_c types for 1D Eulerian meshes
!! \use spatial_discr_rad_m Provides spatial_discr_rad_c type for radial (cylindrical/spherical) discretization
!! \use transport_properties_heterog_m Provides heterogeneous transport properties (dispersion, flux, source terms)

subroutine compute_trans_mat_tpt_trans_Lagr(this)
!> T: transition matrix (tridiagonal, negative semi-definite) [1/T]
!> rows sum = 0 if r=0 (mass conservation constraint)
!> F*dc/dt=T*c+g (semi-discrete transport equation in Lagrangian framework)
    use transport_transient_m, only: transport_1D_transient_c !> Transient 1D transport class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c !> 1D Eulerian mesh types (homogeneous/heterogeneous)
    use spatial_discr_rad_m, only: spatial_discr_rad_c !> Radial spatial discretization class (cylindrical/spherical symmetry)
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c !> Heterogeneous transport properties (dispersion, flux, source terms)
    implicit none !> Enforce explicit variable declaration (no implicit typing)
    
    class(transport_1D_transient_c) :: this !> Transport object containing mesh, properties, and transition matrix to be computed
    
    real(kind=8) :: sign_flux !> Sign of Darcy flux (+1 for positive flow, -1 for negative flow, currently unused in Lagrangian) [-]
    real(kind=8) :: r_i_12 !> Radius at interface i+1/2 for radial symmetry (not currently used) [L]
    real(kind=8) :: a,b !> Auxiliary variables for non-uniform mesh spacing: a = Δxᵢ₋₁/₂, b = Δxᵢ₊₁/₂ [L]
    integer(kind=4) :: i,n !> Loop index i over cells [1..n], total number of cells n [-]
    
    n=this%spatial_discr%Num_targets !> Get total number of spatial nodes (cells or interfaces) from mesh discretization [-]

    if (this%spatial_discr%adapt_ref.eq.1) then !> Check if adaptive mesh refinement flag is enabled (adapt_ref = 1 means mesh has changed) [-]
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super) !> Deallocate old transition matrix arrays (sub, diag, super) to prepare for new mesh
        call this%allocate_trans_mat() !> Reallocate transition matrix arrays with new mesh size (n-1, n, n-1 entries)
    end if
    
    select type (mesh=>this%spatial_discr) !> Polymorphic type selection to determine mesh type (homogeneous/radial) and compute matrix accordingly
    type is (mesh_1D_Euler_homog_c) !> Case: 1D homogeneous Eulerian mesh with uniform spacing Δx (Lagrangian formulation)
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then !> Dimensionless formulation, cell-centered targets (targets_flag=0)
            if (mesh%scheme<3 .and. mesh%scheme>0) then !> Centered finite differences (scheme 1 or 2) for dimensionless Lagrangian case
                this%trans_mat%sub=1d0/(mesh%Delta_x**2) !+ 1d0/(2*mesh%Delta_x) !> Subdiagonal: Tᵢ,ᵢ₋₁ = 1/Δx² (dispersion only, advection commented out) [1/L²]
                this%trans_mat%super=1d0/(mesh%Delta_x**2) !- 1d0/(2*mesh%Delta_x) !> Superdiagonal: Tᵢ,ᵢ₊₁ = 1/Δx² (dispersion only, advection commented out) [1/L²]
            !else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024) scheme - alternative formulation (commented out - not active)
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> First superdiagonal entry
            !    do i=2,n-1 !> Loop over interior cells
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal entries
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal entries
            !    end do
            !    this%trans_mat%sub(n-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Last subdiagonal entry
            else if (mesh%scheme.eq.3) then !> Upwind finite difference scheme (scheme=3) - Lagrangian with dispersion only
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then !> Check if flux has consistent sign (all positive or all negative)
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1)) !> Extract flux sign: +1 if flux > 0, -1 if flux < 0 (not used in Lagrangian) [-]
                    this%trans_mat%sub=1d0/(mesh%Delta_x**2)!+((sign_flux+1d0)/2)*1d0/mesh%Delta_x !> Subdiagonal: T_sub = 1/Δx² (dispersion only, upwind advection commented out) [1/L²]
                    this%trans_mat%super=1d0/(mesh%Delta_x**2)!+((sign_flux-1d0)/2)*1d0/mesh%Delta_x !> Superdiagonal: T_super = 1/Δx² (dispersion only, upwind advection commented out) [1/L²]
                end if
            else !> Invalid scheme number for dimensionless Lagrangian case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported scheme
            end if
        else if (mesh%targets_flag .eq. 0) then !> Dimensional formulation, cell-centered targets (targets_flag=0) in Lagrangian framework
            if (mesh%scheme.eq.1) then !> Traditional centered finite differences (scheme=1) with dimensional dispersion properties
                this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2) !+ this%tpt_props_heterog%flux_cent(2:n)& !> Subdiagonal: Tᵢ,ᵢ₋₁ = Dᵢ/Δx² (dispersion only, advection term commented out) [1/T]
                    !/(2*mesh%Delta_x) !> Advection term not included in Lagrangian formulation (commented out)
                this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) !- & !> Superdiagonal: Tᵢ,ᵢ₊₁ = Dᵢ/Δx² (dispersion only, advection term commented out) [1/T]
                    !this%tpt_props_heterog%flux_cent(1:n-1)/(2*mesh%Delta_x) !> Advection term not included in Lagrangian formulation (commented out)
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024) interface-based scheme (scheme=2) for Lagrangian formulation
                this%trans_mat%super(1)=this%tpt_props_heterog%disp_int(2)/(mesh%Delta_x**2) !- & !> First superdiagonal: D₁₊₁/₂/Δx² using interface dispersion (advection commented out) [1/T]
                    !this%tpt_props_heterog%flux_int(2)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                this%trans_mat_prev%super(1)=this%tpt_props_heterog%disp_int(1)/(mesh%Delta_x**2) !this%trans_mat%super(1) !this%tpt_props_heterog%disp_int(2)/(mesh%Delta_x**2) !+ & !> First subdiagonal: D₁₋₁/₂/Δx² using interface dispersion (advection commented out) [1/T]
                    !this%tpt_props_heterog%flux_int(2)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                do i=2,n-1 !> Loop over interior cells i=2..n-1 to compute matrix entries using interface dispersion values
                    this%trans_mat%sub(i-1)=this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) !+ & !> Subdiagonal: Dᵢ₋₁/₂/Δx² at left interface (advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_int(i)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                    this%trans_mat%super(i)=this%tpt_props_heterog%disp_int(i+1)/(mesh%Delta_x**2) !- & !> Superdiagonal: Dᵢ₊₁/₂/Δx² at right interface (advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_int(i+1)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                    this%trans_mat_prev%sub(i-1)=this%tpt_props_heterog%disp_int(i-1)/(mesh%Delta_x**2) !+ &
                        !this%tpt_props_heterog%disp_int(i-1) !this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) !+ & !> Subdiagonal: Dᵢ₋₁/₂/Δx² at left interface (advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_int(i)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                    this%trans_mat_prev%super(i)=this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) !this%tpt_props_heterog%disp_int(i+1)/(mesh%Delta_x**2) !- & !> Superdiagonal: Dᵢ₊₁/₂/Δx² at right interface (advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_int(i+1)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                end do
                ! this%trans_mat_prev%sub(1)=this%trans_mat_prev%sub(1)*2d0
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%disp_int(n)/(mesh%Delta_x**2) !+ & !> Last subdiagonal: Dₙ₋₁/₂/Δx² using interface dispersion (advection commented out) [1/T]
                    !this%tpt_props_heterog%flux_int(n)/(2*mesh%Delta_x) !> Advection term not included (commented out)
                this%trans_mat_prev%sub(n-1)=this%tpt_props_heterog%disp_int(n-1)/(mesh%Delta_x**2) !this%tpt_props_heterog%disp_int(n)/(mesh%Delta_x**2) !+ & !> Last subdiagonal: Dₙ₋₁/₂/Δx² using interface dispersion (advection commented out) [1/T]
                    !this%tpt_props_heterog%flux_int(n)/(2*mesh%Delta_x) !> Advection term not included (commented out)
            else if (mesh%scheme.eq.3) then !> Upwind scheme (scheme=3) for dimensional Lagrangian case (dispersion only)
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then !> Check if flux has uniform sign across domain (no flow reversal)
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1)) !> Determine flow direction: +1 for forward, -1 for backward (not used in Lagrangian) [-]
                    this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2)!+((sign_flux+1d0)/2)*& !> Subdiagonal: D/Δx² (dispersion only, upwind advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_cent(2:n)/mesh%Delta_x !> Upwind advection not included in Lagrangian (commented out)
                    this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2)!+((sign_flux-1d0)/2)*& !> Superdiagonal: D/Δx² (dispersion only, upwind advection commented out) [1/T]
                        !this%tpt_props_heterog%flux_cent(1:n-1)/mesh%Delta_x !> Upwind advection not included in Lagrangian (commented out)
                end if
            else !> Invalid scheme number for dimensional cell-centered Lagrangian case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported scheme
            end if
        else if (this%dimless.eqv..false. .and. this%tpt_props_heterog%homog_flag.eqv..true.) & !> Edge-centered targets, dimensional formulation, homogeneous properties in Lagrangian framework
            then !> (targets located at cell interfaces rather than centers)
            this%trans_mat%sub=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) !& !> Subdiagonal: Tᵢ,ᵢ₋₁ = D_int,i/Δx² using interface dispersion values (advection commented out) [1/T]
            this%trans_mat_prev%sub(1)=this%trans_mat%sub(1) !> chapuza para primera subdiagonal entry (not used in Lagrangian, set equal to first computed value for consistency)
            this%trans_mat_prev%sub(2:n-1)=this%trans_mat%sub(1:n-2)
            this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) !& !> Superdiagonal: Tᵢ,ᵢ₊₁ = D_int,i/Δx² using interface dispersion values (advection commented out) [1/T]
            this%trans_mat_prev%super(2:n-1)=this%trans_mat%super(1:n-2)
                !- this%tpt_props_heterog%flux_int(2:n)/(2*mesh%Delta_x) !> Advection term not included in Lagrangian formulation (commented out)
        end if
    type is (spatial_discr_rad_c) !> Case: Radial spatial discretization (cylindrical or spherical symmetry) with Lagrangian formulation
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then !> Dimensionless formulation, cell-centered targets for radial Lagrangian (not implemented)
            !if (mesh%scheme.eq.1) then !> Traditional centered finite differences (commented out - not implemented)
            !    this%trans_mat%sub=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal entry (placeholder)
            !    this%trans_mat%super=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal entry (placeholder)
            !else if (mesh%scheme.eq.2) then !> Petchamé & Carrera scheme (commented out - not implemented)
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> First superdiagonal
            !    do i=2,n-1 !> Loop over interior cells
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x) !> Subdiagonal entries
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x) !> Superdiagonal entries
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
            !    error stop "Scheme not implemented yet"
            !end if
        else if ((mesh%targets_flag.eq.0) .and. (this%dimless.eqv..false.)) then !> Dimensional formulation, cell-centered targets for radial Lagrangian geometry
            if (mesh%scheme.eq.1) then !> Traditional centered finite differences (scheme=1) for radial case (commented out - not implemented)
                !this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2) + this%tpt_props_heterog%flux(2:n)/& !> Subdiagonal (placeholder)
                !(2*mesh%Delta_x) !> Continuation
                !this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2) - & !> Superdiagonal (placeholder)
                !this%tpt_props_heterog%flux(1:n-1)/(2*mesh%Delta_x) !> Continuation
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024) interface-based scheme for non-uniform radial mesh with Lagrangian formulation
                    a=5d-1*mesh%Delta_r(1) !> Compute distance to left interface: a = rᵢ - rᵢ₋₁ [L]
                    b=mesh%targets(2)%coord(1)-mesh%targets(1)%coord(1) !> Compute distance to right interface: b = rᵢ₊₁ - rᵢ [L]
                    this%trans_mat%super(1)=2d0*this%tpt_props_heterog%disp_int(2)/(b*(a+b)) !> Superdiagonal: Tᵢ,ᵢ₊₁ = 2Dᵢ₊₁/₂/(b·(a+b)) (dispersion only, no advection) [1/T]
                    !print *, this%trans_mat%super(1) !> Debugging output: display first superdiagonal entry value
                do i=2,n-1 !> Loop over interior radial targets i=2..n-1 to compute Lagrangian transition matrix entries (dispersion only)
                    !a=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i-1)) !> Compute half-spacing to left interface: a = (Δrᵢ + Δrᵢ₋₁)/2 [L]
                    !b=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1)) !> Compute half-spacing to right interface: b = (Δrᵢ + Δrᵢ₊₁)/2 [L]
                    a=mesh%targets(i)%coord(1)-mesh%targets(i-1)%coord(1) !> Compute distance to left interface: a = rᵢ - rᵢ₋₁ [L]
                    !print *, mesh%targets(i+1)%coord(1)
                    b=mesh%targets(i+1)%coord(1)-mesh%targets(i)%coord(1) !> Compute distance to right interface: b = rᵢ₊₁ - rᵢ [L]
                    this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%disp_int(i)/(a*(a+b)) !> Subdiagonal: Tᵢ,ᵢ₋₁ = 2Dᵢ₋₁/₂/(a·(a+b)) (dispersion only, no advection) [1/T]
                    this%trans_mat%super(i)=2d0*this%tpt_props_heterog%disp_int(i+1)/(b*(a+b)) !> Superdiagonal: Tᵢ,ᵢ₊₁ = 2Dᵢ₊₁/₂/(b·(a+b)) (dispersion only, no advection) [1/T]
                    !print *, i, this%trans_mat%sub(i-1), this%trans_mat%super(i) !> Debugging output: display cell index and matrix entries for each interior cell
                end do
                b=5d-1*mesh%Delta_r(n) !> Compute distance to left interface: a = rᵢ - rᵢ₋₁ [L]
                a=mesh%targets(n)%coord(1)-mesh%targets(n-1)%coord(1) !> Compute distance to right interface: b = rᵢ₊₁ - rᵢ [L]
                this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%disp_int(n)/(a*(a+b)) !> Subdiagonal: Tᵢ,ᵢ₋₁ = 2Dᵢ₋₁/₂/(a·(a+b)) (dispersion only, no advection) [1/T]
            else if (mesh%scheme.eq.3) then !> Upwind scheme (scheme=3) for radial Lagrangian case (commented out - not fully implemented)
                !if (minval(this%tpt_props_heterog%flux)>=0d0 .or. maxval(this%tpt_props_heterog%flux)<0d0) then !> Check for uniform flux direction
                !    sign_flux=sign(1d0,this%tpt_props_heterog%flux(1)) !> Extract flux sign
                !    this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*& !> Subdiagonal with upwind (commented out)
                !    this%tpt_props_heterog%flux(2:n)/mesh%Delta_x !> Upwind advection term (commented out)
                !    this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*& !> Superdiagonal with upwind (commented out)
                !    this%tpt_props_heterog%flux(1:n-1)/mesh%Delta_x !> Upwind advection term (commented out)
                !end if
            else !> Invalid scheme number for radial cell-centered Lagrangian case
                error stop "Scheme not implemented yet" !> Terminate program with error message for unsupported radial scheme
            end if
        else if (mesh%targets_flag.eq.1 .and. this%dimless.eqv..false.) then !> Dimensional formulation, edge-centered targets for radial Lagrangian (targets at interfaces)
            this%trans_mat%super(1)=this%tpt_props_heterog%disp_cent(1)/(mesh%Delta_r(1)**2) !> First superdiagonal: T₁,₂ = 2D₁/Δr₁² for edge-centered radial discretization [1/T]
                do i=2,n-1 !> Loop over interior interface nodes i=2..n-1 for edge-centered radial Lagrangian formulation
                    !a=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i-1)) !> Half-spacing computation (commented out - not used in this formulation)
                    !b=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1)) !> Half-spacing computation (commented out - not used in this formulation)
                    this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%disp_cent(i-1)/( & !> Subdiagonal: Tᵢ,ᵢ₋₁ = 2Dᵢ / [Δrᵢ₋₁·(Δrᵢ₋₁+Δrᵢ)] (dispersion only) [1/T]
                        mesh%Delta_r(i-1)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) !> Denominator: product of left cell size and total spacing
                    this%trans_mat%super(i)=2d0*this%tpt_props_heterog%disp_cent(i)/( & !> Superdiagonal: Tᵢ,ᵢ₊₁ = 2Dᵢ / [Δrᵢ·(Δrᵢ₋₁+Δrᵢ)] (dispersion only) [1/T]
                        mesh%Delta_r(i)*(mesh%Delta_r(i-1)+mesh%Delta_r(i))) !> Denominator: product of right cell size and total spacing
                end do
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%disp_cent(n-1)/(mesh%Delta_r(n-1)**2) !> Last subdiagonal: Tₙ,ₙ₋₁ = 2Dₙ/Δrₙ₋₁² for edge-centered radial discretization [1/T]
        else !> Invalid configuration for radial Lagrangian mesh (unsupported target or dimensionless combination)
            error stop "Radial transition matrix not fully implemented yet" !> Terminate program indicating radial Lagrangian implementation is incomplete
        end if
    end select !> End polymorphic mesh type selection
    this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term !> Initialize diagonal with negative source term: Tᵢᵢ = -fᵢ·rᵢ (flag multiplies source) [1/T]
    this%trans_mat%diag(2:n-1)=this%trans_mat%diag(2:n-1)-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1) !> Complete diagonal for interior nodes: Tᵢᵢ = -(Tᵢ,ᵢ₋₁ + Tᵢ,ᵢ₊₁) - r ensuring row sum = -r [1/T]
    this%trans_mat_prev%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term !> Initialize diagonal with negative source term: Tᵢᵢ = -fᵢ·rᵢ (flag multiplies source) [1/T]
    this%trans_mat_prev%diag(2:n-1)=this%trans_mat_prev%diag(2:n-1)-this%trans_mat_prev%sub(1:n-2)-this%trans_mat_prev%super(2:n-1) !> Complete diagonal for interior nodes in previous matrix: Tᵢᵢ = -(Tᵢ,ᵢ₋₁ + Tᵢ,ᵢ₊₁) - r [1/T]
end subroutine !> End of compute_trans_mat_tpt_trans_Lagr subroutine 