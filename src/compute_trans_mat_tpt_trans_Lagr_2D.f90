!> \file compute_trans_mat_tpt_trans_Lagr_2D.f90
!! \brief Compute transition matrix T for transient 2D transport using Lagrangian formulation
!! \details This subroutine constructs the spatial operator (transition matrix) for the transient
!! transport equation in 2D using finite differences with Lagrangian framework (moving with flow).
!!
!! The Lagrangian formulation follows the fluid particles, removing the advection term from the
!! spatial operator. The discretized transport equation has the form:
!! \f[
!!   \mathbf{F} \frac{d\mathbf{c}}{dt} = \mathbf{T} \mathbf{c} + \mathbf{g}
!! \f]
!! where:
!! - \f$\mathbf{F}\f$ is the storage matrix (diagonal with porosity) [-]
!! - \f$\mathbf{T}\f$ is the transition matrix (pentadiagonal, negative semi-definite, symmetric) [1/T]
!! - \f$\mathbf{c}\f$ is the concentration vector [C]
!! - \f$\mathbf{g}\f$ is the external terms vector [C/T]
!!
!! **Key difference from Eulerian formulation:**
!! The Lagrangian approach tracks fluid particles moving with the flow, so the transition matrix
!! contains only dispersion terms (no advection). This makes T symmetric (sub = super for
!! each pair of coupled cells), which enables Cholesky decomposition for matrix inversion.
!!
!! Properties of transition matrix \f$\mathbf{T}\f$:
!! - Pentadiagonal structure (sub2, sub, diag, super, super2) for 2D grid
!! - Symmetric: T_{ij} = T_{ji} (dispersion only, no advection asymmetry)
!! - Negative semi-definite (maximum principle preservation)
!! - Row sums equal zero when source term r = 0 (mass conservation)
!!
!! Matrix bandwidth structure for Nx x Ny grid (row-major ordering):
!! - sub/super: bandwidth 1 (x-direction neighbours within same row)
!! - sub2/super2: bandwidth Nx (y-direction neighbours in adjacent rows)
!!
!! Finite difference schemes:
!! - scheme = 1: Traditional centered finite differences (2nd order, dispersion only)
!!
!! \param[in,out] this Transport object (transport_2D_transient_c)

subroutine compute_trans_mat_tpt_trans_Lagr_2D(this)
!> T: transition matrix (pentadiagonal, negative semi-definite, symmetric) [1/T]
!> rows sum = 0 if r=0 (mass conservation constraint)
!> F*dc/dt=T*c+g (semi-discrete transport equation in Lagrangian framework)
    use transport_transient_m, only: transport_2D_transient_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c
    use transport_properties_heterog_m, only: tpt_props_heterog_2D_c
    implicit none

    class(transport_2D_transient_c) :: this

    integer(kind=4) :: i,j,n,Num_cells

    n=this%spatial_discr%Num_targets
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag

    if (this%spatial_discr%adapt_ref.eq.1) then
        deallocate(this%trans_mat%sub, this%trans_mat%diag, this%trans_mat%super, &
                   this%trans_mat%sub2, this%trans_mat%super2)
        call this%allocate_trans_mat()
    end if

    select type (mesh=>this%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then
            !> Dimensionless Lagrangian: dispersion only (no advection)
            if (mesh%scheme<3 .and. mesh%scheme>0) then
                this%trans_mat%sub=1d0/(mesh%Delta_x**2)
                this%trans_mat%super=1d0/(mesh%Delta_x**2)
            else if (mesh%scheme.eq.3) then
                this%trans_mat%sub=1d0/(mesh%Delta_x**2)
                this%trans_mat%super=1d0/(mesh%Delta_x**2)
            else
                error stop "Scheme not implemented yet"
            end if
        else if (mesh%targets_flag .eq. 0) then
            !> Dimensional Lagrangian, cell-centered targets
            if (mesh%scheme.eq.1) then
                !> sub/super (x-direction): size Ny*(Nx-1), flat index = (i-1)*(Nx-1)+j
                !> Lagrangian: dispersion only (D_xx/Δx²), advection removed
                do i=1,mesh%Num_cells_y
                    do j=1,mesh%Num_cells_x-1
                        !> sub flat index: (i-1)*(Nx-1)+j couples cell (i,j+1) to (i,j)
                        !> Dispersion: D_xx at x-interface (j+1) in row i, divided by Δx² [1/T]
                        this%trans_mat%sub((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j+1,i,1)/(mesh%Delta_x**2)
                    end do
                end do
                do i=1,mesh%Num_cells_y
                    do j=1,mesh%Num_cells_x-1
                        !> super flat index: (i-1)*(Nx-1)+j couples cell (i,j) to (i,j+1)
                        !> Dispersion: D_xx at x-interface (j+1) in row i, divided by Δx² [1/T]
                        !> Note: sub = super for the same interface (symmetric matrix)
                        this%trans_mat%super((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j+1,i,1)/(mesh%Delta_x**2)
                    end do
                end do
                !> sub2/super2 (y-direction): size (Ny-1)*Nx, flat index = (i-1)*Nx+j
                !> Lagrangian: dispersion only (D_yy/Δy²), advection removed
                do i=1,mesh%Num_cells_y-1
                    do j=1,mesh%Num_cells_x
                        !> sub2: couples cell (i+1,j) to cell (i,j) above
                        !> Dispersion: D_yy at y-interface (i+1) in column j, divided by Δy² [1/T]
                        this%trans_mat%sub2((j-1)*(mesh%Num_cells_y-1)+i)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i+1,2)/(mesh%Delta_y**2)
                    end do
                end do
                do i=1,mesh%Num_cells_y-1
                    do j=1,mesh%Num_cells_x
                        !> super2: couples cell (i,j) to cell (i+1,j) below
                        !> Dispersion: D_yy at y-interface (i+1) in column j, divided by Δy² [1/T]
                        !> Note: sub2 = super2 for the same interface (symmetric matrix)
                        this%trans_mat%super2((j-1)*(mesh%Num_cells_y-1)+i)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i+1,2)/(mesh%Delta_y**2)
                    end do
                end do
            else
                error stop "Scheme not implemented yet for 2D Lagrangian"
            end if
        else if (.not.this%dimless .and. this%tpt_props_heterog%homog_flag) then
            !> Edge-centered targets, dimensional, homogeneous properties
            !> Lagrangian: dispersion only (D/Δx², D/Δy²), no advection
            this%trans_mat%sub=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2)
            this%trans_mat%super=this%tpt_props_heterog%disp_int(2:n)/(mesh%Delta_x**2)
        end if
        !> Diagonal: initialise with negative source term, then enforce row-sum = -r
        this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term
        !> --- interior nodes: rows 2..Ny-1, cols 2..Nx-1 ---
        do i=2,mesh%Num_cells_y-1
            do j=2,mesh%Num_cells_x-1
                this%trans_mat%diag((i-1)*mesh%Num_cells_x+j)= &
                    this%trans_mat%diag((i-1)*mesh%Num_cells_x+j) &
                    - this%trans_mat%sub( (i-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                    - this%trans_mat%super( (i-1)*(mesh%Num_cells_x-1)+j ) &
                    - this%trans_mat%sub2( (j-1)*(mesh%Num_cells_y-1)+(i-1) ) &
                    - this%trans_mat%super2( (j-1)*(mesh%Num_cells_y-1)+i )
            end do
        end do
        !> --- top boundary row (i=1): cols 2..Nx-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat%diag(j)= &
                this%trans_mat%diag(j) &
                - this%trans_mat%sub( j-1 ) &
                - this%trans_mat%super( j ) &
                - this%trans_mat%super2( (j-1)*(mesh%Num_cells_y-1)+1 )
        end do
        !> --- bottom boundary row (i=Ny): cols 2..Nx-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j)= &
                this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j) &
                - this%trans_mat%sub( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                - this%trans_mat%super( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+j ) &
                - this%trans_mat%sub2( (j-1)*(mesh%Num_cells_y-1)+(mesh%Num_cells_y-1) )
        end do
        !> --- left boundary column (j=1): rows 2..Ny-1 ---
        do i=2,mesh%Num_cells_y-1
            this%trans_mat%diag((i-1)*mesh%Num_cells_x+1)= &
                this%trans_mat%diag((i-1)*mesh%Num_cells_x+1) &
                - this%trans_mat%super( (i-1)*(mesh%Num_cells_x-1)+1 ) &
                - this%trans_mat%sub2( i-1 ) &
                - this%trans_mat%super2( i )
        end do
        !> --- right boundary column (j=Nx): rows 2..Ny-1 ---
        do i=2,mesh%Num_cells_y-1
            this%trans_mat%diag(i*mesh%Num_cells_x)= &
                this%trans_mat%diag(i*mesh%Num_cells_x) &
                - this%trans_mat%sub( (i-1)*(mesh%Num_cells_x-1)+(mesh%Num_cells_x-1) ) &
                - this%trans_mat%sub2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+(i-1) ) &
                - this%trans_mat%super2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+i )
        end do
        !> --- corner (1,1): top-left ---
        this%trans_mat%diag(1)= &
            this%trans_mat%diag(1) &
            - this%trans_mat%super( 1 ) &
            - this%trans_mat%super2( 1 )
        !> --- corner (1,Nx): top-right ---
        this%trans_mat%diag(mesh%Num_cells_x)= &
            this%trans_mat%diag(mesh%Num_cells_x) &
            - this%trans_mat%sub( mesh%Num_cells_x-1 ) &
            - this%trans_mat%super2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+1 )
        !> --- corner (Ny,1): bottom-left ---
        this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+1)= &
            this%trans_mat%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+1) &
            - this%trans_mat%super( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+1 ) &
            - this%trans_mat%sub2( mesh%Num_cells_y-1 )
        !> --- corner (Ny,Nx): bottom-right ---
        this%trans_mat%diag(mesh%Num_cells_y*mesh%Num_cells_x)= &
            this%trans_mat%diag(mesh%Num_cells_y*mesh%Num_cells_x) &
            - this%trans_mat%sub( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+(mesh%Num_cells_x-1) ) &
            - this%trans_mat%sub2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+(mesh%Num_cells_y-1) )

    end select
    !> Debug: verify rows of trans_mat sum to 0
    block
        integer(kind=4) :: k, col_k, row_k, Nx_dbg, Ny_dbg
        real(kind=8) :: rsum
        select type (mesh=>this%spatial_discr)
        type is (mesh_2D_Euler_homog_c)
            Nx_dbg = mesh%Num_cells_x
            Ny_dbg = mesh%Num_cells_y
        end select
        do k = 1, n
            col_k = mod(k-1, Nx_dbg) + 1
            row_k = (k-1)/Nx_dbg + 1
            rsum = this%trans_mat%diag(k)
            if (col_k > 1)      rsum = rsum + this%trans_mat%sub((row_k-1)*(Nx_dbg-1) + col_k - 1)
            if (col_k < Nx_dbg) rsum = rsum + this%trans_mat%super((row_k-1)*(Nx_dbg-1) + col_k)
            if (row_k > 1)      rsum = rsum + this%trans_mat%sub2((col_k-1)*(Ny_dbg-1) + row_k - 1)
            if (row_k < Ny_dbg) rsum = rsum + this%trans_mat%super2((col_k-1)*(Ny_dbg-1) + row_k)
            !print '(A,I4,A,ES14.6)', '[DEBUG Lagr 2D] trans_mat row ', k, ' sum = ', rsum
            if (abs(rsum) > 1d-6) then
                error stop "trans_mat row does not sum to 0 in Lagr 2D"
            end if
        end do
    end block
    !> Compute trans_mat_prev for isochrone mesh (shifted interface indices by -1)
    !> In the isochrone Lagrangian framework, particles shift by one cell per time step,
    !> so the previous mesh uses interface dispersion shifted by -1 in each direction:
    !>   x-direction: disp_tensor_int(j+1,i,1) -> disp_tensor_int(j,i,1)
    !>   y-direction: disp_tensor_int(j,i+1,2) -> disp_tensor_int(j,i,2)
    select type (mesh=>this%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        if (mesh%targets_flag .eq. 0 .and. .not.this%dimless) then
            if (mesh%scheme.eq.1) then
                !> sub_prev (x-direction): shifted index j instead of j+1
                do i=1,mesh%Num_cells_y
                    do j=1,mesh%Num_cells_x-1
                        this%trans_mat_prev%sub((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i,1)/(mesh%Delta_x**2)
                    end do
                end do
                !> super_prev (x-direction): shifted index j instead of j+1 (symmetric)
                do i=1,mesh%Num_cells_y
                    do j=1,mesh%Num_cells_x-1
                        this%trans_mat_prev%super((i-1)*(mesh%Num_cells_x-1)+j)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i,1)/(mesh%Delta_x**2)
                    end do
                end do
                !> sub2_prev (y-direction): shifted index i instead of i+1
                do i=1,mesh%Num_cells_y-1
                    do j=1,mesh%Num_cells_x
                        this%trans_mat_prev%sub2((j-1)*(mesh%Num_cells_y-1)+i)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i,2)/(mesh%Delta_y**2)
                    end do
                end do
                !> super2_prev (y-direction): shifted index i instead of i+1 (symmetric)
                do i=1,mesh%Num_cells_y-1
                    do j=1,mesh%Num_cells_x
                        this%trans_mat_prev%super2((j-1)*(mesh%Num_cells_y-1)+i)=&
                            this%tpt_props_heterog%disp_tensor_int(&
                            j,i,2)/(mesh%Delta_y**2)
                    end do
                end do
            end if
        else
            !> Non-isochrone cases: copy trans_mat into trans_mat_prev
            this%trans_mat_prev%sub=this%trans_mat%sub
            this%trans_mat_prev%super=this%trans_mat%super
            this%trans_mat_prev%sub2=this%trans_mat%sub2
            this%trans_mat_prev%super2=this%trans_mat%super2
        end if
        !> Diagonal of trans_mat_prev: row-sum = -r (mass conservation)
        this%trans_mat_prev%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term
        !> --- interior nodes: rows 2..Ny-1, cols 2..Nx-1 ---
        do i=2,mesh%Num_cells_y-1
            do j=2,mesh%Num_cells_x-1
                this%trans_mat_prev%diag((i-1)*mesh%Num_cells_x+j)= &
                    this%trans_mat_prev%diag((i-1)*mesh%Num_cells_x+j) &
                    - this%trans_mat_prev%sub( (i-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                    - this%trans_mat_prev%super( (i-1)*(mesh%Num_cells_x-1)+j ) &
                    - this%trans_mat_prev%sub2( (j-1)*(mesh%Num_cells_y-1)+(i-1) ) &
                    - this%trans_mat_prev%super2( (j-1)*(mesh%Num_cells_y-1)+i )
            end do
        end do
        !> --- top boundary row (i=1): cols 2..Nx-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat_prev%diag(j)= &
                this%trans_mat_prev%diag(j) &
                - this%trans_mat_prev%sub( j-1 ) &
                - this%trans_mat_prev%super( j ) &
                - this%trans_mat_prev%super2( (j-1)*(mesh%Num_cells_y-1)+1 )
        end do
        !> --- bottom boundary row (i=Ny): cols 2..Nx-1 ---
        do j=2,mesh%Num_cells_x-1
            this%trans_mat_prev%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j)= &
                this%trans_mat_prev%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+j) &
                - this%trans_mat_prev%sub( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+(j-1) ) &
                - this%trans_mat_prev%super( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+j ) &
                - this%trans_mat_prev%sub2( (j-1)*(mesh%Num_cells_y-1)+(mesh%Num_cells_y-1) )
        end do
        !> --- left boundary column (j=1): rows 2..Ny-1 ---
        do i=2,mesh%Num_cells_y-1
            this%trans_mat_prev%diag((i-1)*mesh%Num_cells_x+1)= &
                this%trans_mat_prev%diag((i-1)*mesh%Num_cells_x+1) &
                - this%trans_mat_prev%super( (i-1)*(mesh%Num_cells_x-1)+1 ) &
                - this%trans_mat_prev%sub2( i-1 ) &
                - this%trans_mat_prev%super2( i )
        end do
        !> --- right boundary column (j=Nx): rows 2..Ny-1 ---
        do i=2,mesh%Num_cells_y-1
            this%trans_mat_prev%diag(i*mesh%Num_cells_x)= &
                this%trans_mat_prev%diag(i*mesh%Num_cells_x) &
                - this%trans_mat_prev%sub( (i-1)*(mesh%Num_cells_x-1)+(mesh%Num_cells_x-1) ) &
                - this%trans_mat_prev%sub2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+(i-1) ) &
                - this%trans_mat_prev%super2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+i )
        end do
        !> --- corner (1,1): top-left ---
        this%trans_mat_prev%diag(1)= &
            this%trans_mat_prev%diag(1) &
            - this%trans_mat_prev%super( 1 ) &
            - this%trans_mat_prev%super2( 1 )
        !> --- corner (1,Nx): top-right ---
        this%trans_mat_prev%diag(mesh%Num_cells_x)= &
            this%trans_mat_prev%diag(mesh%Num_cells_x) &
            - this%trans_mat_prev%sub( mesh%Num_cells_x-1 ) &
            - this%trans_mat_prev%super2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+1 )
        !> --- corner (Ny,1): bottom-left ---
        this%trans_mat_prev%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+1)= &
            this%trans_mat_prev%diag((mesh%Num_cells_y-1)*mesh%Num_cells_x+1) &
            - this%trans_mat_prev%super( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+1 ) &
            - this%trans_mat_prev%sub2( mesh%Num_cells_y-1 )
        !> --- corner (Ny,Nx): bottom-right ---
        this%trans_mat_prev%diag(mesh%Num_cells_y*mesh%Num_cells_x)= &
            this%trans_mat_prev%diag(mesh%Num_cells_y*mesh%Num_cells_x) &
            - this%trans_mat_prev%sub( (mesh%Num_cells_y-1)*(mesh%Num_cells_x-1)+(mesh%Num_cells_x-1) ) &
            - this%trans_mat_prev%sub2( (mesh%Num_cells_x-1)*(mesh%Num_cells_y-1)+(mesh%Num_cells_y-1) )
    end select
end subroutine
