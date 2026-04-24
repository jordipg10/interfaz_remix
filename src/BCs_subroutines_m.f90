!> \file BCs_subroutines_m.f90
!> \brief Module for imposing boundary conditions on PDE transition matrices and source terms
!> \details
!>   Provides subroutines to modify the transition matrix T, source term vector g,
!>   and boundary matrix bd_mat for various BC combinations:
!>   - Dirichlet–Dirichlet
!>   - Dirichlet–Neumann
!>   - Neumann–Neumann (homogeneous)
!>   - Robin–Neumann (homogeneous)
!>   - Robin–Dirichlet
!>
!>   Each subroutine dispatches on PDE type (transport_1D, transport_2D, flow)
!>   and mesh type (homogeneous Euler, radial, etc.) to apply the correct stencil
!>   modifications at the inflow (index 1) and outflow (index n) boundaries.
module BCs_subroutines_m
    use PDE_m, only: PDE_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c
    use transport_m, only: transport_1D_c
    use transport_transient_m, only: transport_1D_transient_c, transport_2D_transient_c
    use diffusion_transient_m, only: diffusion_1D_transient_c
    use flow_transient_m, only: flow_transient_c
    implicit none
    save
    contains

    !> \brief Imposes Dirichlet boundary conditions at both inflow and outflow
    !> \details Modifies T, bd_mat, and source_term_PDE for Dirichlet BCs
    !> \param[in,out] this PDE object (polymorphic)
    subroutine Dirichlet_BCs_PDE(this)
        implicit none
        class(PDE_c) :: this

        integer(kind=4) :: n
        real(kind=8) :: q_inf
        real(kind=8) :: dx2_inv, D1_dx2

        n = this%spatial_discr%Num_targets

        select type (this)
        type is (transport_1D_transient_c)
            q_inf = this%BCs%flux_inf
            select type (mesh => this%spatial_discr)
            type is (mesh_1D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                if (mesh%targets_flag == 1) then
                    !> Edge-centred finite differences
                    this%trans_mat%diag(1)  = 0d0
                    this%trans_mat%super(1) = 0d0
                    this%trans_mat%diag(n)  = 0d0
                    this%trans_mat%sub(n-1) = 0d0
                    this%trans_mat_prev%diag(1)  = 0d0
                    this%trans_mat_prev%super(1) = 0d0
                    D1_dx2 = this%tpt_props_heterog%disp_cent(n-1) * dx2_inv
                    this%trans_mat_prev%diag(n) = this%trans_mat_prev%diag(n) &
                        - (this%tpt_props_heterog%disp_cent(n-1) + this%tpt_props_heterog%disp_cent(n-2)) * dx2_inv
                    this%bd_mat_prev(2) = D1_dx2
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 1) then
                    !> Cell-centred FD, traditional scheme
                    call apply_Dirichlet_CCFD_trad_1D(this, mesh, n, dx2_inv)
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 2) then
                    !> Cell-centred FD, Petchamé & Carrera (2024)
                    call apply_Dirichlet_CCFD_prop_1D(this, mesh, n, dx2_inv)
                else
                    error stop "BCs not implemented yet"
                end if
            type is (spatial_discr_rad_c)
                if (.not. this%dimless) then
                    if (mesh%targets_flag == 0 .and. mesh%scheme == 2) then
                        call apply_Dirichlet_radial_1D(this, mesh, n, q_inf)
                    end if
                end if
            end select

        type is (transport_2D_transient_c)
            select type (mesh => this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                if (mesh%targets_flag == 1) then
                    !> Edge-centred FD (2D)
                    this%trans_mat%diag(1)  = 0d0
                    this%trans_mat%super(1) = 0d0
                    this%trans_mat%diag(n)  = 0d0
                    this%trans_mat%sub(n-1) = 0d0
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 1) then
                    !> CCFD traditional (2D)
                    call apply_Dirichlet_CCFD_trad_2D(this, mesh, n, dx2_inv)
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 2) then
                    !> CCFD proposed (2D)
                    call apply_Dirichlet_CCFD_prop_2D(this, mesh, n, dx2_inv)
                else
                    error stop "BCs not implemented yet"
                end if
            end select

        type is (flow_transient_c)
            select type (mesh => this%spatial_discr)
            type is (spatial_discr_rad_c)
                if (this%dimless) then
                    call apply_Dirichlet_flow_radial(this, mesh, n)
                end if
            end select
        end select
    end subroutine


    !> \brief Imposes Dirichlet inflow & Neumann outflow boundary conditions
    !> \param[in,out] this PDE object (polymorphic)
    subroutine Dirichlet_Neumann_BCs_PDE(this)
        implicit none
        class(PDE_c) :: this

        integer(kind=4) :: n
        real(kind=8) :: q_inf, dx2_inv

        n = this%spatial_discr%Num_targets

        select type (this)
        type is (transport_1D_transient_c)
            q_inf = this%BCs%flux_inf
            select type (mesh => this%spatial_discr)
            type is (mesh_1D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                if (mesh%targets_flag == 1) then
                    !> ECFD
                    this%trans_mat%diag(1)  = 0d0
                    this%trans_mat%super(1) = 0d0
                    this%trans_mat%sub(n-1) = 2d0 * this%tpt_props_heterog%disp_cent(n) * dx2_inv
                    this%trans_mat%diag(n)  = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 1) then
                    !> CCFD traditional
                    call apply_Dir_Neu_CCFD_trad_1D(this, mesh, n, dx2_inv)
                else if (mesh%targets_flag == 0 .and. mesh%scheme == 2) then
                    !> CCFD proposed
                    call apply_Dir_Neu_CCFD_prop_1D(this, mesh, n, dx2_inv, q_inf, &
                        this%tpt_props_heterog%disp_int(1))
                else
                    error stop "BCs not implemented yet"
                end if
            type is (spatial_discr_rad_c)
                if (.not. this%dimless) then
                    if (mesh%targets_flag == 0 .and. mesh%scheme == 2) then
                        call apply_Dir_Neu_radial_1D(this, mesh, n, q_inf)
                    else if (mesh%targets_flag == 1) then
                        this%trans_mat%diag(1)  = 0d0
                        this%trans_mat%super(1) = 0d0
                        this%bd_mat(1) = 0d0
                        this%trans_mat%sub(n-1) = 2d0 * this%tpt_props_heterog%disp_int(n) / (mesh%Delta_r(n-1)**2)
                        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
                    end if
                end if
            end select
        end select
    end subroutine


    !> \brief Imposes Neumann homogeneous BCs at both inflow and outflow
    !> \param[in,out] this PDE object (polymorphic)
    subroutine Neumann_homog_BCs(this)
        implicit none
        class(PDE_c) :: this

        integer(kind=4) :: n

        n = this%spatial_discr%Num_targets

        select type (this)
        type is (transport_1D_transient_c)
            if (this%spatial_discr%targets_flag == 1) then
                !> ECFD: double the boundary dispersion entries
                this%trans_mat%super(1) = 2d0 * this%tpt_props_heterog%disp_int(1) &
                    / (this%spatial_discr%get_cell_size(1)**2)
                this%trans_mat%sub(n-1) = 2d0 * this%tpt_props_heterog%disp_int(n) &
                    / (this%spatial_discr%get_cell_size(n-1)**2)
            else if (this%spatial_discr%targets_flag /= 0) then
                error stop "BCs not implemented yet"
            end if
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
            this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
        end select
    end subroutine


    !> \brief Imposes Robin BC at inflow & Neumann homogeneous BC at outflow
    !> \param[in,out] this PDE object (polymorphic)
    subroutine Robin_Neumann_homog_BCs(this)
        implicit none
        class(PDE_c) :: this

        integer(kind=4) :: n
        real(kind=8) :: q_cent_1, q_cent_n, q_inf, q_out, D_1, D_n, q_int_1
        real(kind=8) :: a, b, c, d, dx2_inv

        n = this%spatial_discr%Num_targets

        select type (this)
        type is (transport_1D_transient_c)
            q_cent_1 = this%tpt_props_heterog%flux_cent(1)
            q_int_1  = this%tpt_props_heterog%flux_int(1)
            q_cent_n = this%tpt_props_heterog%flux_cent(n - this%spatial_discr%targets_flag)
            q_inf    = this%BCs%flux_inf
            q_out    = this%BCs%flux_out
            D_1      = this%tpt_props_heterog%disp_int(1)
            D_n      = this%tpt_props_heterog%disp_int(n)

            select type (mesh => this%spatial_discr)
            type is (mesh_1D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                a = 2d0*D_1 - q_inf*mesh%Delta_x
                b = q_inf*mesh%Delta_x + 2d0*D_1
                c = q_cent_1*mesh%Delta_x + 2d0*D_1
                d = 2d0*D_1 - q_cent_1*mesh%Delta_x

                if (mesh%scheme == 1 .and. mesh%targets_flag == 0) then
                    call apply_Robin_Neu_CCFD_trad_1D(this, mesh, n, dx2_inv, a, b, c, d, q_inf)
                else if (mesh%scheme == 2 .and. mesh%targets_flag == 0) then
                    call apply_Robin_Neu_CCFD_prop_1D(this, mesh, n, dx2_inv, q_inf, D_1)
                else if (mesh%targets_flag == 1) then
                    call apply_Robin_Neu_ECFD_1D(this, mesh, n, dx2_inv, q_inf, q_int_1, D_1, D_n)
                end if
            type is (spatial_discr_rad_c)
                if (.not. this%dimless) then
                    if (mesh%scheme == 2 .and. mesh%targets_flag == 0 &
                        .and. (.not. this%tpt_props_heterog%homog_flag)) then
                        a = mesh%Delta_r(1)
                        b = mesh%Delta_r(1) + mesh%Delta_r(2)
                        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                            - q_inf/a + this%tpt_props_heterog%flux_int(2)/b &
                            - 4d0*this%tpt_props_heterog%disp_int(2)/(b**2)
                        this%trans_mat%super(1) = -this%tpt_props_heterog%flux_int(2)/b &
                            + 4d0*this%tpt_props_heterog%disp_int(2)/(b**2)
                        this%bd_mat(1) = q_inf / a
                    end if
                end if
            end select
            this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)

        type is (transport_2D_transient_c)
            q_cent_1 = this%tpt_props_heterog%flux_cent(1)
            q_int_1  = this%tpt_props_heterog%flux_int(1)
            q_cent_n = this%tpt_props_heterog%flux_cent(n - this%spatial_discr%targets_flag)
            q_inf    = this%BCs%flux_inf
            q_out    = this%BCs%flux_out
            D_1      = this%tpt_props_heterog%disp_int(1)
            D_n      = this%tpt_props_heterog%disp_int(n)

            select type (mesh => this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                a = 2d0*D_1 - q_inf*mesh%Delta_x
                b = q_inf*mesh%Delta_x + 2d0*D_1
                c = q_cent_1*mesh%Delta_x + 2d0*D_1
                d = 2d0*D_1 - q_cent_1*mesh%Delta_x

                if (mesh%scheme == 1 .and. mesh%targets_flag == 0) then
                    call apply_Robin_Neu_CCFD_trad_2D(this, mesh, n, dx2_inv, a, b, c, d, q_inf, D_1, D_n)
                else if (mesh%scheme == 2 .and. mesh%targets_flag == 0) then
                    call apply_Robin_Neu_CCFD_prop_2D(this, mesh, n, dx2_inv, q_inf, D_1, D_n)
                else if (mesh%targets_flag == 1) then
                    call apply_Robin_Neu_ECFD_2D(this, mesh, n, dx2_inv, q_inf, q_int_1, D_1, D_n)
                end if
            end select
        end select
    end subroutine


    !> \brief Imposes Robin BC at inflow & Dirichlet BC at outflow
    !> \param[in,out] this PDE object (polymorphic)
    subroutine Robin_Dirichlet_BCs(this)
        implicit none
        class(PDE_c) :: this

        integer(kind=4) :: n
        real(kind=8) :: q_cent_1, q_cent_n, q_inf, q_out, c_inf, dx2_inv

        n = this%spatial_discr%Num_targets

        select type (this)
        type is (transport_1D_transient_c)
            q_cent_1 = this%tpt_props_heterog%flux_cent(1)
            q_cent_n = this%tpt_props_heterog%flux_cent(n)
            q_inf    = this%BCs%flux_inf
            q_out    = this%BCs%flux_out
            c_inf    = this%BCs%conc_inf

            select type (mesh => this%spatial_discr)
            type is (mesh_1D_Euler_homog_c)
                dx2_inv = 1d0 / (mesh%Delta_x**2)
                if (this%tpt_props_heterog%homog_flag .and. mesh%targets_flag == 0) then
                    call apply_Robin_Dir_CCFD_homog_1D(this, mesh, n, dx2_inv, q_cent_1, q_cent_n, q_inf, q_out)
                else if (mesh%targets_flag == 1 .and. this%tpt_props_heterog%homog_flag) then
                    call apply_Robin_Dir_ECFD_homog_1D(this, mesh, n, dx2_inv, q_inf)
                end if
            end select

        end select
    end subroutine


    !=========================================================================
    ! Private helper subroutines
    !=========================================================================

    !> \brief Apply Dirichlet BCs for 1D CCFD traditional scheme
    subroutine apply_Dirichlet_CCFD_trad_1D(this, mesh, n, dx2_inv)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv
        real(kind=8) :: half_dx_inv, D1, Dn

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        D1 = this%tpt_props_heterog%disp_cent(1)
        Dn = this%tpt_props_heterog%disp_cent(n)

        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            - this%tpt_props_heterog%flux_cent(1) * half_dx_inv - 2d0*D1*dx2_inv
        this%trans_mat%super(1) = -this%tpt_props_heterog%flux_cent(1)*half_dx_inv + D1*dx2_inv
        this%trans_mat%diag(n) = this%trans_mat%diag(n) &
            + this%tpt_props_heterog%flux_cent(n)*half_dx_inv - 3d0*Dn*dx2_inv
        this%trans_mat%sub(n-1) = this%tpt_props_heterog%flux_cent(n)*half_dx_inv + Dn*dx2_inv

        this%bd_mat(1) = (this%tpt_props_heterog%flux_cent(1) + 2d0*D1/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        this%bd_mat(2) = -(this%tpt_props_heterog%flux_cent(n) - 2d0*Dn/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
    end subroutine


    !> \brief Apply Dirichlet BCs for 1D CCFD proposed scheme (Petchamé & Carrera 2024)
    subroutine apply_Dirichlet_CCFD_prop_1D(this, mesh, n, dx2_inv)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv
        real(kind=8) :: half_dx_inv, D1, D2, Dn, Dn1

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        D1  = this%tpt_props_heterog%disp_int(1)
        D2  = this%tpt_props_heterog%disp_int(2)
        Dn  = this%tpt_props_heterog%disp_int(n)
        Dn1 = this%tpt_props_heterog%disp_int(n+1)

        if (.not. this%Lagr_flag) then
            !> Eulerian
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                - this%tpt_props_heterog%flux_int(1)/mesh%Delta_x &
                + this%tpt_props_heterog%flux_int(2)*half_dx_inv &
                - 2d0*D1*dx2_inv - D2*dx2_inv
            this%trans_mat_prev%diag(1) = this%trans_mat%diag(1)
            this%trans_mat%super(1) = -this%tpt_props_heterog%flux_int(2)*half_dx_inv + D2*dx2_inv
            this%trans_mat_prev%super(1) = this%trans_mat%super(1)

            this%trans_mat%diag(n) = this%trans_mat%diag(n) &
                + this%tpt_props_heterog%flux_int(n+1)/mesh%Delta_x &
                - this%tpt_props_heterog%flux_int(n)*half_dx_inv &
                - 2d0*Dn1*dx2_inv - Dn*dx2_inv
            this%trans_mat_prev%diag(n) = this%trans_mat%diag(n)
            this%trans_mat%sub(n-1) = this%tpt_props_heterog%flux_int(n)*half_dx_inv + Dn*dx2_inv
            this%trans_mat_prev%sub(n-1) = this%trans_mat%sub(n-1)

            this%bd_mat(1) = (this%tpt_props_heterog%flux_int(1) + 2d0*D1/mesh%Delta_x) / mesh%Delta_x
            this%bd_mat_prev(1) = this%bd_mat(1)
            this%bd_mat(2) = (-this%tpt_props_heterog%flux_int(n+1) + 2d0*Dn1/mesh%Delta_x) / mesh%Delta_x
            this%bd_mat_prev(2) = this%bd_mat(2)
        else
            !> Lagrangian
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - 2d0*D1*dx2_inv - D2*dx2_inv
            this%trans_mat%super(1) = D2 * dx2_inv
            this%trans_mat_prev%super(1) = 2d0*D1*dx2_inv
            this%trans_mat_prev%diag(1) = this%trans_mat_prev%diag(1) - this%trans_mat_prev%super(1)
            this%trans_mat_prev%diag(2) = this%trans_mat_prev%diag(2) - this%trans_mat_prev%sub(1)
            this%trans_mat_prev%sub(1) = this%trans_mat_prev%sub(1) * 2d0
            this%bd_mat(1) = 2d0*D1*dx2_inv

            this%trans_mat%diag(n) = this%trans_mat%diag(n) - 2d0*Dn1*dx2_inv - Dn*dx2_inv
            this%trans_mat_prev%diag(n) = this%trans_mat_prev%diag(n) &
                - this%tpt_props_heterog%disp_int(n)*dx2_inv &
                - this%tpt_props_heterog%disp_int(n-1)*dx2_inv
            this%trans_mat%sub(n-1) = Dn * dx2_inv
            this%trans_mat_prev%sub(n-1) = this%tpt_props_heterog%disp_int(n-1) * dx2_inv

            this%bd_mat(2) = 2d0*Dn1*dx2_inv
            this%bd_mat_prev(2) = Dn * dx2_inv
        end if
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
    end subroutine


    !> \brief Apply Dirichlet BCs for 1D radial mesh (scheme 2)
    subroutine apply_Dirichlet_radial_1D(this, mesh, n, q_inf)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(spatial_discr_rad_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: q_inf
        real(kind=8) :: dr1, dr_sum, D1_coeff, D2_coeff

        dr1 = mesh%Delta_r(1)
        dr_sum = 3d0*dr1 + mesh%Delta_r(2)
        D1_coeff = 8d0*this%tpt_props_heterog%disp_int(1) / (dr1 * (&
            2d0*mesh%Delta_r(1) + mesh%Delta_r(2)))
        D2_coeff = 8d0*this%tpt_props_heterog%disp_int(2) / (dr_sum * (&
            dr1 + mesh%Delta_r(2)))

        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            - q_inf/dr1 + this%tpt_props_heterog%flux_int(2)/(dr1 + mesh%Delta_r(2)) &
            - D1_coeff - D2_coeff
        this%trans_mat%super(1) = -this%tpt_props_heterog%flux_int(2)/(dr1 + mesh%Delta_r(2)) + D2_coeff
        this%bd_mat(1) = q_inf/dr1 + D1_coeff
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
    end subroutine


    !> \brief Apply Dirichlet BCs for 2D CCFD traditional scheme
    subroutine apply_Dirichlet_CCFD_trad_2D(this, mesh, n, dx2_inv)
        implicit none
        type(transport_2D_transient_c), intent(inout) :: this
        type(mesh_2D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv
        real(kind=8) :: half_dx_inv

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)

        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            - this%tpt_props_heterog%flux_cent(1)*half_dx_inv - 3d0*this%tpt_props_heterog%disp_int(1)*dx2_inv
        this%trans_mat%super(1) = -this%tpt_props_heterog%flux_cent(1)*half_dx_inv &
            + this%tpt_props_heterog%disp_int(1)*dx2_inv
        this%trans_mat%diag(n) = this%trans_mat%diag(n) &
            + this%tpt_props_heterog%flux_cent(n)*half_dx_inv - 3d0*this%tpt_props_heterog%disp_int(n)*dx2_inv
        this%trans_mat%sub(n-1) = this%tpt_props_heterog%flux_cent(n)*half_dx_inv &
            + this%tpt_props_heterog%disp_int(n)*dx2_inv

        this%bd_mat(1) = (this%tpt_props_heterog%flux_cent(1) &
            + 2d0*this%tpt_props_heterog%disp_int(1)/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        this%bd_mat(2) = -(this%tpt_props_heterog%flux_cent(n) &
            - 2d0*this%tpt_props_heterog%disp_int(n)/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
    end subroutine


    !> \brief Apply Dirichlet BCs for 2D CCFD proposed scheme
    subroutine apply_Dirichlet_CCFD_prop_2D(this, mesh, n, dx2_inv)
        implicit none
        type(transport_2D_transient_c), intent(inout) :: this
        type(mesh_2D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv
        real(kind=8) :: half_dx_inv, D1, Dn1

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        D1  = this%tpt_props_heterog%disp_int(1)
        Dn1 = this%tpt_props_heterog%disp_int(n+1)

        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            - this%tpt_props_heterog%flux_int(1)/mesh%Delta_x &
            + this%tpt_props_heterog%flux_int(2)*half_dx_inv &
            - 3d0*D1*dx2_inv
        this%trans_mat%diag(n) = this%trans_mat%diag(n) &
            + this%tpt_props_heterog%flux_int(n+1)/mesh%Delta_x &
            - this%tpt_props_heterog%flux_int(n)*half_dx_inv &
            - 3d0*Dn1*dx2_inv

        this%bd_mat(1) = (this%tpt_props_heterog%flux_int(1) + 2d0*D1/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        this%bd_mat(2) = -(this%tpt_props_heterog%flux_int(n+1) - 2d0*Dn1/mesh%Delta_x) / mesh%Delta_x
        this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
    end subroutine


    !> \brief Apply Dirichlet BCs for dimensionless radial flow
    subroutine apply_Dirichlet_flow_radial(this, mesh, n)
        implicit none
        type(flow_transient_c), intent(inout) :: this
        type(spatial_discr_rad_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        integer(kind=4) :: dim_m1
        real(kind=8) :: r_min_D_pow, r_max_D_pow

        dim_m1 = mesh%dim - 1
        r_min_D_pow = mesh%r_min_D**dim_m1
        r_max_D_pow = mesh%r_max_D**dim_m1

        !> Inflow
        this%trans_mat%diag(1) = -4d0 * (r_min_D_pow/mesh%Delta_r_D(1) &
            + (mesh%r_min_D + mesh%Delta_r_D(1))**dim_m1 / sum(mesh%Delta_r_D(1:2)))
        this%bd_mat(1) = 4d0 * r_min_D_pow / mesh%Delta_r_D(1)
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%bd_mat(1) * this%BCs%head_inf_D
        !> Outflow
        this%trans_mat%diag(n) = -4d0 * (r_max_D_pow/mesh%Delta_r_D(n) &
            + (mesh%r_max_D - mesh%Delta_r_D(n))**dim_m1 / sum(mesh%Delta_r_D(n-1:n)))
        this%bd_mat(2) = 4d0 * r_max_D_pow / mesh%Delta_r_D(n)
        this%source_term_PDE(n) = this%source_term_PDE(n) + this%bd_mat(2) * this%BCs%head_out_D
    end subroutine


    !> \brief Apply Dirichlet-Neumann BCs for 1D CCFD traditional scheme (no inflow flux parameters)
    subroutine apply_Dir_Neu_CCFD_trad_1D(this, mesh, n, dx2_inv)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv

        if (this%Lagr_flag) then
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
        end if
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
    end subroutine


    !> \brief Apply Dirichlet-Neumann BCs for 1D CCFD proposed scheme
    subroutine apply_Dir_Neu_CCFD_prop_1D(this, mesh, n, dx2_inv, q_inf, D_1)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf, D_1
        real(kind=8) :: half_dx_inv, Dn

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        Dn = this%tpt_props_heterog%disp_int(n)

        if (.not. this%Lagr_flag) then
            !> Eulerian
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                - q_inf/mesh%Delta_x + this%tpt_props_heterog%flux_int(2)*half_dx_inv &
                - this%tpt_props_heterog%disp_int(2)*dx2_inv
            this%bd_mat(1) = q_inf / mesh%Delta_x
        else
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
            this%trans_mat_prev%super(1) = 0d0
            this%trans_mat_prev%diag(2) = this%trans_mat_prev%diag(2) + this%trans_mat_prev%sub(1)
            this%trans_mat_prev%sub(1) = 0d0
            this%bd_mat_prev(2) = Dn * dx2_inv
            this%trans_mat_prev%diag(n) = this%trans_mat_prev%diag(n) &
                - this%trans_mat_prev%sub(n-1) - this%bd_mat_prev(2)
        end if
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 1D CCFD traditional scheme
    subroutine apply_Robin_Neu_CCFD_trad_1D(this, mesh, n, dx2_inv, a, b, c, d, q_inf)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, a, b, c, d, q_inf

        if (.not. this%Lagr_flag) then
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                + a*c / (b * 2d0 * mesh%Delta_x**2) &
                - 2d0*this%tpt_props_heterog%disp_int(1)*dx2_inv
            this%trans_mat%super(1) = d / (2d0 * mesh%Delta_x**2)
            this%bd_mat(1) = q_inf * c / (b * mesh%Delta_x)
        else
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
        end if
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 1D CCFD proposed scheme
    subroutine apply_Robin_Neu_CCFD_prop_1D(this, mesh, n, dx2_inv, q_inf, D_1)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf, D_1
        real(kind=8) :: half_dx_inv, Dn

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        Dn = this%tpt_props_heterog%disp_int(n)

        if (this%Lagr_flag .eqv. .false.) then
            !> Eulerian
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                - q_inf/mesh%Delta_x + this%tpt_props_heterog%flux_int(2)*half_dx_inv &
                - this%tpt_props_heterog%disp_int(2)*dx2_inv
            this%trans_mat_prev%diag(1) = this%trans_mat%diag(1)
            this%bd_mat(1) = q_inf / mesh%Delta_x
            this%bd_mat_prev(1) = this%bd_mat(1)
            this%trans_mat_prev%super(1) = this%trans_mat%super(1)
        else
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
            this%trans_mat_prev%super(1) = 0d0
            this%trans_mat_prev%diag(1) = this%trans_mat_prev%diag(1) - this%trans_mat_prev%super(1)
            this%trans_mat_prev%diag(2) = this%trans_mat_prev%diag(2) + this%trans_mat_prev%sub(1)
            this%trans_mat_prev%sub(1) = 0d0
            this%bd_mat_prev(2) = Dn * dx2_inv
            this%trans_mat_prev%diag(n) = this%trans_mat_prev%diag(n) &
                - this%trans_mat_prev%sub(n-1) - this%bd_mat_prev(2)
        end if
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
        if (.not. this%Lagr_flag) then
            this%trans_mat_prev%diag(n) = this%trans_mat%diag(n)
            this%trans_mat_prev%sub(n-1) = this%trans_mat%sub(n-1)
        end if
    end subroutine


    !> \brief Apply Dirichlet-Neumann BCs for 1D radial mesh (scheme 2)
    subroutine apply_Dir_Neu_radial_1D(this, mesh, n, q_inf)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(spatial_discr_rad_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: q_inf
        real(kind=8) :: dr1, dr_sum, D1_coeff, D2_coeff

        dr1 = mesh%Delta_r(1)
        dr_sum = 2d0*dr1 + mesh%Delta_r(2)
        D1_coeff = 8d0*this%tpt_props_heterog%disp_int(1) / (dr1 * (&
            2d0*mesh%Delta_r(1) + mesh%Delta_r(2)))
        D2_coeff = 8d0*this%tpt_props_heterog%disp_int(2) / (dr_sum * (&
            dr1 + mesh%Delta_r(2)))

        if (.not. this%Lagr_flag) then
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                - q_inf/dr1 + this%tpt_props_heterog%flux_int(2)/(dr1 + mesh%Delta_r(2)) &
                - D1_coeff - D2_coeff
            this%trans_mat%super(1) = -this%tpt_props_heterog%flux_int(2)/(dr1 + mesh%Delta_r(2)) + D2_coeff
            this%bd_mat(1) = q_inf / dr1
            this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        else
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - D1_coeff - D2_coeff
            this%trans_mat%super(1) = D2_coeff
            this%bd_mat(1) = D1_coeff
            this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        end if
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 2D CCFD traditional scheme
    subroutine apply_Robin_Neu_CCFD_trad_2D(this, mesh, n, dx2_inv, a, b, c, d, q_inf, D_1, D_n)

        implicit none
        type(transport_2D_transient_c), intent(inout) :: this
        type(mesh_2D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, a, b, c, d, q_inf, D_1, D_n
        integer(kind=4) :: i, j
        real(kind=8) :: dy2_inv

        dy2_inv = 1d0 / (mesh%Delta_y**2)

        if (.not. this%Lagr_flag) then
            do i = 1, mesh%Num_cells_y
                this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) = &
                    this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) &
                    + a*c/(b*2d0*mesh%Delta_x**2) - 2d0*D_1*dx2_inv
                this%trans_mat%super((i-1)*(mesh%Num_cells_x-1) + 1) = d / (2d0*mesh%Delta_x**2)
            end do
            this%bd_mat(1) = q_inf*c / (b*mesh%Delta_x)
        end if

        ! !> Apply row-sum correction at x-boundaries for each y-row
        ! do j = 1, mesh%Num_cells_y
        !     this%trans_mat%super((j-1)*(mesh%Num_cells_x-1) + 1) = &
        !         this%tpt_props_heterog%disp_tensor_int(2, j, 1) * dx2_inv
        !     this%trans_mat%diag((j-1)*mesh%Num_cells_x + 1) = &
        !         this%trans_mat%diag((j-1)*mesh%Num_cells_x + 1) &
        !         - this%trans_mat%super((j-1)*(mesh%Num_cells_x-1) + 1)
        !     this%trans_mat%sub(j*(mesh%Num_cells_x-1)) = &
        !         this%tpt_props_heterog%disp_tensor_int(mesh%Num_cells_x, j, 1) * dx2_inv
        !     this%trans_mat%diag(j*mesh%Num_cells_x) = &
        !         this%trans_mat%diag(j*mesh%Num_cells_x) &
        !         - this%trans_mat%sub(j*(mesh%Num_cells_x-1))
        ! end do

        ! !> Apply y-direction diagonal corrections for left/right boundary columns at interior y-rows
        ! do i = 2, mesh%Num_cells_y - 1
        !     !> Left column (j=1): subtract sub2 (upper) and super2 (lower) from diagonal
        !     this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) = &
        !         this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) &
        !         - this%trans_mat%sub2((i-2)*mesh%Num_cells_x + 1) &
        !         - this%trans_mat%super2((i-1)*mesh%Num_cells_x + 1)
        !     !> Right column (j=Nx): subtract sub2 (upper) and super2 (lower) from diagonal
        !     this%trans_mat%diag(i*mesh%Num_cells_x) = &
        !         this%trans_mat%diag(i*mesh%Num_cells_x) &
        !         - this%trans_mat%sub2((i-2)*mesh%Num_cells_x + mesh%Num_cells_x) &
        !         - this%trans_mat%super2((i-1)*mesh%Num_cells_x + mesh%Num_cells_x)
        ! end do

        ! !> Apply row-sum correction at y-boundaries for each x-column
        ! do i = 1, mesh%Num_cells_x
        !     this%trans_mat%super2(i) = this%tpt_props_heterog%disp_tensor_int(i, 2, 2) * dy2_inv
        !     this%trans_mat%diag(i) = this%trans_mat%diag(i) - this%trans_mat%super2(i)
        !     this%trans_mat%sub2((mesh%Num_cells_y-2)*mesh%Num_cells_x + i) = &
        !         this%tpt_props_heterog%disp_tensor_int(i, mesh%Num_cells_y, 2) * dy2_inv
        !     this%trans_mat%diag(n - mesh%Num_cells_x + i) = &
        !         this%trans_mat%diag(n - mesh%Num_cells_x + i) &
        !         - this%trans_mat%sub2((mesh%Num_cells_y-2)*mesh%Num_cells_x + i)
        ! end do
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 2D CCFD proposed scheme
    subroutine apply_Robin_Neu_CCFD_prop_2D(this, mesh, n, dx2_inv, q_inf, D_1, D_n)
        implicit none
        type(transport_2D_transient_c), intent(inout) :: this
        type(mesh_2D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf, D_1, D_n
        real(kind=8) :: half_dx_inv

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)

        if (.not. this%Lagr_flag) then
            this%trans_mat%diag(1) = this%trans_mat%diag(1) &
                - q_inf/mesh%Delta_x + this%tpt_props_heterog%flux_int(2)*half_dx_inv &
                - this%tpt_props_heterog%disp_int(2)*dx2_inv
            this%bd_mat(1) = q_inf / mesh%Delta_x
        else
            this%trans_mat%diag(1) = this%trans_mat%diag(1) - this%trans_mat%super(1)
            this%trans_mat_prev%super(1) = 0d0
            this%trans_mat_prev%diag(2) = this%trans_mat_prev%diag(2) + this%trans_mat_prev%sub(1)
            this%trans_mat_prev%sub(1) = 0d0
            this%bd_mat_prev(2) = D_n * dx2_inv
            this%trans_mat_prev%diag(n) = this%trans_mat_prev%diag(n) &
                - this%trans_mat_prev%sub(n-1) - this%bd_mat_prev(2)
        end if

        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
    end subroutine


    !> \brief Apply Robin-Dirichlet BCs for 1D CCFD homogeneous properties
    subroutine apply_Robin_Dir_CCFD_homog_1D(this, mesh, n, dx2_inv, q_cent_1, q_cent_n, q_inf, q_out)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_cent_1, q_cent_n, q_inf, q_out
        real(kind=8) :: half_dx_inv, D1, Dn, D1_int

        half_dx_inv = 1d0 / (2d0 * mesh%Delta_x)
        D1     = this%tpt_props_heterog%disp_cent(1)
        Dn     = this%tpt_props_heterog%disp_cent(n)
        D1_int = this%tpt_props_heterog%disp_int(1)

        this%trans_mat%super(1) = -q_cent_1*half_dx_inv + D1*dx2_inv
        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            + (2d0*D1_int - q_inf*mesh%Delta_x)*(q_cent_1*mesh%Delta_x + 2d0*D1) &
            / (2d0*q_inf*mesh%Delta_x**3 + 4d0*D1_int*mesh%Delta_x**2) &
            - 2d0*D1*dx2_inv
        this%bd_mat(1) = q_inf*(q_cent_1*mesh%Delta_x + 2d0*D1) &
            / (q_inf*mesh%Delta_x**2 + 2d0*D1_int*mesh%Delta_x)
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)

        if (mesh%scheme == 1) then
            !> Traditional scheme at outflow
            this%trans_mat%diag(n) = this%trans_mat%diag(n) &
                + q_cent_n*half_dx_inv - 3d0*Dn*dx2_inv
            this%trans_mat%sub(n-1) = q_cent_n*half_dx_inv + Dn*dx2_inv
            this%bd_mat(2) = -(q_cent_n - 2d0*Dn/mesh%Delta_x) / mesh%Delta_x
            this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
        else if (mesh%scheme == 2) then
            !> Proposed scheme at outflow
            this%trans_mat%diag(n) = this%trans_mat%diag(n) &
                + q_out/mesh%Delta_x &
                - this%tpt_props_heterog%flux_int(n)*half_dx_inv &
                - 2d0*this%tpt_props_heterog%disp_int(n+1)*dx2_inv &
                - this%tpt_props_heterog%disp_int(n)*dx2_inv
            this%trans_mat%sub(n-1) = this%tpt_props_heterog%flux_int(n)*half_dx_inv &
                + this%tpt_props_heterog%disp_int(n)*dx2_inv
            this%bd_mat(2) = -(q_out - 2d0*this%tpt_props_heterog%disp_int(n+1)/mesh%Delta_x) / mesh%Delta_x
            this%source_term_PDE(n) = this%source_term_PDE(n) + this%BCs%conc_out * this%bd_mat(2)
        end if
    end subroutine


    !> \brief Apply Robin-Dirichlet BCs for 1D ECFD homogeneous properties
    subroutine apply_Robin_Dir_ECFD_homog_1D(this, mesh, n, dx2_inv, q_inf)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf
        real(kind=8) :: D1_int, two_D1_dx2

        D1_int     = this%tpt_props_heterog%disp_int(1)
        two_D1_dx2 = 2d0 * D1_int * dx2_inv

        this%trans_mat%super(1) = two_D1_dx2
        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            - q_inf*(this%tpt_props_heterog%flux_int(1)*mesh%Delta_x + 2d0*D1_int)/(D1_int*mesh%Delta_x) &
            - two_D1_dx2
        this%bd_mat(1) = q_inf*(this%tpt_props_heterog%flux_int(1)*mesh%Delta_x + 2d0*D1_int) &
            / (D1_int*mesh%Delta_x)
        this%source_term_PDE(1) = this%source_term_PDE(1) + this%BCs%conc_inf * this%bd_mat(1)
        this%trans_mat%diag(n) = 0d0
        this%trans_mat%sub(n-1) = 0d0
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 1D ECFD (edge-centred finite differences)
    subroutine apply_Robin_Neu_ECFD_1D(this, mesh, n, dx2_inv, q_inf, q_int_1, D_1, D_n)
        implicit none
        type(transport_1D_transient_c), intent(inout) :: this
        type(mesh_1D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf, q_int_1, D_1, D_n
        real(kind=8) :: a, b, c, d, two_D1_dx2

        !> Robin inflow: same a,b,c,d coefficients as CCFD-trad but using interface dispersion
        a = 2d0*D_1 - q_inf*mesh%Delta_x
        b = q_inf*mesh%Delta_x + 2d0*D_1
        c = q_int_1*mesh%Delta_x + 2d0*D_1
        d = 2d0*D_1 - q_int_1*mesh%Delta_x

        two_D1_dx2 = 2d0 * D_1 * dx2_inv

        this%trans_mat%super(1) = d / (2d0 * mesh%Delta_x**2)
        this%trans_mat%diag(1) = this%trans_mat%diag(1) &
            + a*c / (b * 2d0 * mesh%Delta_x**2) &
            - two_D1_dx2
        this%bd_mat(1) = q_inf * c / (b * mesh%Delta_x)

        !> Neumann homogeneous outflow: double the sub-diagonal (image-cell reflection)
        this%trans_mat%sub(n-1) = 2d0 * D_n * dx2_inv
        this%trans_mat%diag(n) = this%trans_mat%diag(n) - this%trans_mat%sub(n-1)
    end subroutine


    !> \brief Apply Robin-Neumann BCs for 2D ECFD (edge-centred finite differences)
    subroutine apply_Robin_Neu_ECFD_2D(this, mesh, n, dx2_inv, q_inf, q_int_1, D_1, D_n)
        implicit none
        type(transport_2D_transient_c), intent(inout) :: this
        type(mesh_2D_Euler_homog_c), intent(in) :: mesh
        integer(kind=4), intent(in) :: n
        real(kind=8), intent(in) :: dx2_inv, q_inf, q_int_1, D_1, D_n
        real(kind=8) :: a, b, c, d
        integer(kind=4) :: i

        !> Robin inflow: same a,b,c,d coefficients as CCFD-trad but using interface dispersion
        a = 2d0*D_1 - q_inf*mesh%Delta_x
        b = q_inf*mesh%Delta_x + 2d0*D_1
        c = q_int_1*mesh%Delta_x + 2d0*D_1
        d = 2d0*D_1 - q_int_1*mesh%Delta_x

        !> Apply Robin inflow at x=0 boundary for each y-row
        do i = 1, mesh%Num_cells_y
            this%trans_mat%super((i-1)*(mesh%Num_cells_x-1) + 1) = d / (2d0 * mesh%Delta_x**2)
            this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) = &
                this%trans_mat%diag((i-1)*mesh%Num_cells_x + 1) &
                + a*c / (b * 2d0 * mesh%Delta_x**2) &
                - 2d0*D_1*dx2_inv
        end do
        this%bd_mat(1) = q_inf * c / (b * mesh%Delta_x)

        !> Neumann homogeneous outflow: double sub-diagonal at x=L boundary for each y-row
        do i = 1, mesh%Num_cells_y
            this%trans_mat%sub(i*(mesh%Num_cells_x-1)) = 2d0 * D_n * dx2_inv
            this%trans_mat%diag(i*mesh%Num_cells_x) = &
                this%trans_mat%diag(i*mesh%Num_cells_x) &
                - this%trans_mat%sub(i*(mesh%Num_cells_x-1))
        end do
    end subroutine

end module
