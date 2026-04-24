!> \file transport_stab_params_m.f90
!> \brief Transport (advection-dispersion) stability parameter module
!> \details 
!>   Defines stability parameters for advection-dispersion-reaction (ADE) transport.
!>   Extends diffusion stability parameters with advection constraints.
!>
!>   Governing Equation (1D Transport):
!>   \f[
!>     \phi \frac{\partial c}{\partial t} = \frac{\partial}{\partial x}\left(D \frac{\partial c}{\partial x}\right) - \frac{\partial(qc)}{\partial x} + r
!>   \f]
!>
!>   Stability Criteria:
!>   
!>   1. **Courant Number** (Advection Stability):
!>      \f[
!>        Co = \frac{|q| \Delta t}{\phi \Delta x} \leq 1
!>      \f]
!>      - Limits advective transport to one cell per time step
!>      - Prevents numerical oscillations in explicit schemes
!>
!>   2. **Peclet Number** (Spatial Resolution):
!>      \f[
!>        Pe = \frac{\Delta x}{\alpha_L} = \frac{|q| \Delta x}{D} \leq 2
!>      \f]
!>      - Ensures mesh resolves dispersion length scale
!>      - Controls numerical dispersion vs. physical dispersion
!>
!>   3. **Fourier Number** (Diffusion Stability, inherited):
!>      \f[
!>        Fo = \frac{D \Delta t}{\phi \Delta x^2} \leq 0.5
!>      \f]
!>
!>   Class Hierarchy:
!>   - stab_params_tpt_1D_c extends stab_params_diff_c
!>   - Adds Courant and Peclet to inherited Fourier number
!>
!>   Usage:
!>   \code{.f90}
!>   type(stab_params_tpt_1D_c) :: stab_params
!>   call stab_params%compute_stab_params(props, mesh, dt)
!>   if (stab_params%Courant > 1.0) then
!>       ! Reduce time step or refine mesh
!>   end if
!>   \endcode
!>
!>   \see stab_params_diff_c, transport_properties_heterog_m, spatial_discr_1D_m
module transport_stab_params_m
    use diff_stab_params_m, only: stab_params_diff_c
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c, tpt_props_heterog_2D_c
    use spatial_discr_m, only: spatial_discr_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c
    use properties_m, only: props_c
    implicit none
    save
    private
    !> \class stab_params_tpt_1D_c
    !> \brief 1D transport stability parameters class
    !> \extends stab_params_diff_c
    !> \details 
    !>   Encapsulates stability constraints for advection-dispersion transport.
    !>   Inherits diffusion stability (Fourier number, critical time step).
    !>   Adds advection-specific constraints (Courant, Peclet).
    type, public, extends(stab_params_diff_c) :: stab_params_tpt_1D_c
        !> \var Courant
        !> Courant number: Co = |q|Δt/(φΔx) [-]
        !> Advection stability constraint (should be ≤ 1 for explicit schemes)
        real(kind=8) :: Courant
        
        !> \var Peclet
        !> Cell Peclet number: Pe = Δx/α_L [-]
        !> Spatial resolution constraint (should be ≤ 2 to avoid numerical dispersion)
        real(kind=8) :: cell_Peclet
    contains
        !> Compute all transport stability parameters (Fourier, Courant, Peclet)
        procedure :: compute_stab_params_tpt_1D
        procedure :: compute_cell_Peclet_1D
    end type

    type, public, extends(stab_params_tpt_1D_c) :: stab_params_tpt_2D_c
        !> Additional 2D-specific stability parameters can be added here
        real(kind=8) :: Courant_y  !< Courant number in y-direction
        !real(kind=8) :: Peclet_y   !< Peclet number in y-direction
    contains
        !> Compute all transport stability parameters in 2D (Fourier, Courant, Peclet)
        procedure :: compute_stab_params_tpt_2D
        procedure :: compute_cell_Peclet_2D
    end type
    contains
        !> \brief Compute transport stability parameters (Fourier, Courant, Peclet)
        !> \param[in,out] this      Stability parameters object to update
        !> \param[in]     props_obj Transport properties object
        !> \param[in]     mesh      Spatial discretization object
        !> \param[in]     time_step Time step size [T]
        !> \details 
        !>   Computes all stability constraints for advection-dispersion transport.
        !>   Inherits Fourier number computation, adds Courant and Peclet numbers.
        !>
        !>   Computed Parameters:
        !>
        !>   1. **Fourier Number** (from parent class):
        !>      \f[
        !>        Fo = \frac{D \Delta t}{\phi \Delta x^2}
        !>      \f]
        !>
        !>   2. **Courant Number** (advection stability):
        !>      \f[
        !>        Co = \frac{|q| \Delta t}{\phi \Delta x}
        !>      \f]
        !>      Violation (Co > 1): Prints warning, continues simulation
        !>
        !>   3. **Peclet Number** (spatial resolution):
        !>      \f[
        !>        Pe = \frac{\Delta x}{\alpha_L}
        !>      \f]
        !>      Violation (Pe > 2): Prints warning, continues simulation
        !>
        !>   Homogeneous Case:
        !>   - Uses single representative values (first cell)
        !>   - Computes parameters once (uniform properties)
        !>
        !>   Heterogeneous Case:
        !>   - Loops over all cells to find maximum values
        !>   - Co_max = max(|q_i| Δt / (φ_i Δx_i))
        !>   - Pe_max = max(Δx_i / α_L)
        !>   - Conservative approach (worst-case constraint)
        !>
        !>   Assumptions:
        !>   - Mesh and time discretization can be uniform or heterogeneous
        !>   - Transport properties pre-computed (flux, dispersion, porosity)
        !>
        !>   \warning Stability violations print warnings but do not halt execution
        subroutine compute_stab_params_tpt_1D(this,props_obj,mesh,time_step)
            implicit none
            class(stab_params_tpt_1D_c) :: this                    !< Stability parameters object (modified)
            class(tpt_props_heterog_1D_c), intent(in) :: props_obj             !< Properties object (transport)
            class(spatial_discr_c), intent(in) :: mesh          !< Spatial discretization
            real(kind=8), intent(in) :: time_step               !< Time step Δt [T]
            
            real(kind=8) :: D,phi,q                             !< Dispersion [L²/T], porosity [-], flux [L/T]
            real(kind=8) :: Courant,Peclet,Courant_max,Peclet_max,beta_max,beta  !< Stability numbers (local and maximum)
            real(kind=8), parameter :: epsilon=1d-12            !< Tolerance for floating-point comparisons [-]
            integer(kind=4) :: i                                !< Loop index
            
            !call this%stab_params_diff_c%compute_stab_params_diff_1D(props_obj,mesh,time_step)  !< Compute diffusion stability (Fourier, Δt_crit)
            
            select type (mesh)                             !< Polymorphic dispatch on properties type
            type is (mesh_1D_Euler_homog_c)                       !< Case: transport properties
                !call are_tpt_props_homog(props_obj)             !< Check spatial homogeneity
                ! if (props_obj%homog_flag.eqv..true.) then       !< Case: homogeneous properties
                !     phi=props_obj%porosity(1)                   !< Use porosity from first cell
                !     D=props_obj%disp_cent(1)              !< Use dispersion from first cell center
                !     q=props_obj%flux_cent(1)                    !< Use flux from first cell center
                !     this%Delta_t_crit=phi*mesh%get_cell_size()**2/(2d0*D)  !< Critical time step: Δt_crit = φΔx²/(2D)
                !     this%Courant=abs(q)*time_step/(phi*mesh%get_cell_size())  !< Courant number: Co = |q|Δt/(φΔx)
                !     if (this%Courant>1d0) then                  !< Check Courant condition
                !         print *, this%Courant                   !< Print violation value
                !         print *, "Courant condition violated"   !< Print warning message
                !         !error stop  "Courant condition violated"  !< Optional: halt execution
                !     end if
                !     this%Peclet=abs(q)*mesh%get_max_cell_size()/D  !< Peclet number: Pe = |q|Δx/D = Δx/α_L
                !     if (this%Peclet>2d0) then                   !< Check Peclet condition
                !         print *, this%Peclet                    !< Print violation value
                !         print *, "Peclet condition violated"    !< Print warning message
                !         !error stop  "Peclet condition violated"  !< Optional: halt execution
                !     end if
                !else                                            !< Case: heterogeneous properties
                    Courant_max=props_obj%flux_cent(1)*time_step/(props_obj%porosity(1)*mesh%Delta_x)  !< Initialize with first cell
                    beta_max=props_obj%disp_cent(1)*time_step/(props_obj%porosity(1)*mesh%Delta_x**2)  !< Initialize Peclet with first cell value: Pe = Δx/α_L
                    do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                        Courant=props_obj%flux_cent(i)*time_step/(props_obj%porosity(i)*mesh%Delta_x)  !< Compute local Courant
                        beta=props_obj%disp_cent(i)*time_step/(props_obj%porosity(i)*mesh%Delta_x**2)  !< Compute local Peclet: Pe = Δx/α_L
                        if (Courant>Courant_max) then           !< Check if current exceeds maximum
                            Courant_max=Courant                 !< Update maximum
                        end if
                        if (beta>beta_max) then                 !< Check if current exceeds maximum
                            beta_max=beta                       !< Update maximum
                        end if
                    end do
                    this%Courant=Courant_max                    !< Store maximum Courant number
                    this%beta=beta_max                     !< Store maximum Peclet number
                    if (this%Courant>1d0) then                  !< Check Courant condition
                        write(*,'(A,ES10.3,A)') ' WARNING: Courant number Co = ', this%Courant, ' > 1.0'
                        write(*,'(A)') '   Reduce Delta_t or refine mesh (decrease Delta_x)'
                    end if
                    if (this%beta>0.5d0+1d-10) then               !< Check dispersion condition (with tolerance)
                        write(*,'(A,ES10.3,A)') ' WARNING: Dispersion number beta = ', this%beta, ' > 0.5'
                        write(*,'(A)') '   Reduce Delta_t or refine mesh (decrease Delta_x)'
                    end if
                    this%cell_Peclet=mesh%Delta_x/props_obj%long_dispersivity  !< Initialize with first cell: Pe = Δx/α_L
                    ! do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                    !     Peclet=mesh%Delta_x/props_obj%long_dispersivity  !< Compute local Peclet
                    !     if (Peclet>Peclet_max) then             !< Check if current exceeds maximum
                    !         Peclet_max=Peclet                   !< Update maximum
                    !     end if
                    ! end do
                    ! this%Peclet=Peclet_max                      !< Store maximum Peclet number
                    if (this%cell_Peclet>2d0) then                   !< Check Peclet condition
                        write(*,'(A,ES10.3,A)') ' WARNING: Cell Peclet number Pe = ', this%cell_Peclet, ' > 2.0'
                        write(*,'(A)') '   Refine mesh (decrease Delta_x) or increase dispersivity'
                    end if
                !end if
            type is (mesh_1D_Euler_heterog_c)                       !< Case: transport properties
                Courant_max=props_obj%flux_cent(1)*time_step/(props_obj%porosity(1)*mesh%Delta_x(1))  !< Initialize with first cell
                beta_max=props_obj%disp_cent(1)*time_step/(props_obj%porosity(1)*mesh%Delta_x(1)**2)  !< Initialize with first cell
                do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                    Courant=props_obj%flux_cent(i)*time_step/(props_obj%porosity(i)*mesh%Delta_x(i))  !< Compute local Courant
                    beta=props_obj%disp_cent(i)*time_step/(props_obj%porosity(i)*mesh%Delta_x(i)**2)  !< Compute local Peclet: Pe = Δx/α_L
                    if (Courant>Courant_max) then           !< Check if current exceeds maximum
                        Courant_max=Courant                 !< Update maximum
                    end if
                    if (beta>beta_max) then                 !< Check if current exceeds maximum
                        beta_max=beta                       !< Update maximum
                    end if
                end do
                this%beta=beta_max                     !< Store maximum beta number
                this%Courant=Courant_max                    !< Store maximum Courant number
                if (this%Courant>1d0) then                  !< Check Courant condition
                    write(*,'(A,ES10.3,A)') ' WARNING: Courant number Co = ', this%Courant, ' > 1.0'
                    write(*,'(A)') '   Reduce Delta_t or refine mesh (decrease Delta_x)'
                end if
                this%cell_Peclet=maxval(mesh%Delta_x)/props_obj%long_dispersivity  !< Initialize with first cell: Pe = Δx/α_L
                ! do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                !     Peclet=mesh%Delta_x(i)/props_obj%long_dispersivity  !< Compute local Peclet
                !     if (Peclet>Peclet_max) then             !< Check if current exceeds maximum
                !         Peclet_max=Peclet                   !< Update maximum 
                !     end if
                ! end do
                if (this%cell_Peclet>2d0) then                   !< Check Peclet condition
                    write(*,'(A,ES10.3,A)') ' WARNING: Cell Peclet number Pe = ', this%cell_Peclet, ' > 2.0'
                    write(*,'(A)') '   Refine mesh (decrease Delta_x) or increase dispersivity'
                end if
            end select            
        end subroutine

        subroutine compute_stab_params_tpt_2D(this,props_obj,mesh,time_step)
            implicit none
            class(stab_params_tpt_2D_c) :: this                    !< Stability parameters object (modified)
            class(tpt_props_heterog_2D_c), intent(in) :: props_obj             !< Properties object (transport)
            class(spatial_discr_c), intent(in) :: mesh          !< Spatial discretization
            real(kind=8), intent(in) :: time_step               !< Time step Δt [T]
            
            real(kind=8) :: Courant,Peclet,Courant_max,beta,beta_max          !< Stability numbers (local and maximum)
            integer(kind=4) :: i                                !< Loop index
            !call this%compute_stab_params_diff_2D(props_obj,mesh,time_step)
            !call this%compute_stab_params_tpt_1D(props_obj,mesh,time_step)
            
            !> Additional computations for Courant_y and Peclet_y would go here
            select type (mesh)
            type is (mesh_2D_Euler_homog_c)
                !> Homogeneous case computations for y-direction
                Courant_max=props_obj%flux_cent_y(1)*time_step/(&
                    props_obj%porosity(1)*mesh%Delta_y)  !< Initialize with first cell
                beta_max=props_obj%disp_cent(1)*time_step*mesh%sq_hypot/(&
                    props_obj%porosity(1)*mesh%get_max_cell_size()**2)  !< Initialize Peclet with first cell value: Pe = Δx/α_L
                do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                    Courant=props_obj%flux_cent_y(i)*time_step/(&
                        props_obj%porosity(i)*mesh%Delta_y)  !< Compute local Courant
                    beta=props_obj%disp_cent(i)*time_step*mesh%sq_hypot/(&
                        props_obj%porosity(i)*mesh%get_max_cell_size()**2)
                    if (Courant>Courant_max) then           !< Check if current exceeds maximum
                        Courant_max=Courant                 !< Update maximum
                    end if
                    if (beta>beta_max) then                 !< Check if current exceeds maximum
                        beta_max=beta                       !< Update maximum
                    end if
                end do
                this%beta=beta_max                     !< Store maximum beta number
                this%Courant_y=Courant_max                    !< Store maximum Courant number for y-direction
                if (this%Courant_y>1d0-this%Courant) then                  !< Check Courant condition
                    write(*,'(A,ES10.3,A,ES10.3,A)') ' WARNING: Courant_y = ', this%Courant_y, &
                        ' > 1 - Courant_x = ', 1d0-this%Courant, ' (Co_x + Co_y > 1)'
                    write(*,'(A)') '   Reduce Delta_t or refine mesh in y-direction'
                end if
                if (this%beta>0.5d0+1d-10) then               !< Check dispersion condition (with tolerance)
                    write(*,'(A,ES10.3,A)') ' WARNING: Dispersion number beta = ', this%beta, ' > 0.5'
                    write(*,'(A)') '   Reduce Delta_t or refine mesh'
                end if
                this%cell_Peclet = this%compute_cell_Peclet_2D(mesh%Delta_x,mesh%Delta_y,props_obj%long_dispersivity)
                if (this%cell_Peclet>2d0) then                   !< Check Peclet condition
                    write(*,'(A,ES10.3,A)') ' WARNING: Cell Peclet number Pe = ', this%cell_Peclet, ' > 2.0'
                    write(*,'(A)') '   Refine mesh or increase dispersivity'
                end if
                ! this%Peclet_y=mesh%Delta_y/props_obj%long_dispersivity  !< Initialize with first cell: Pe = Δx/α_L
                ! if (this%Peclet_y>2d0-this%Peclet) then                   !< Check Peclet condition
                !     print *, this%Peclet_y                    !< Print violation value
                !     print *, "Peclet condition in y-direction violated"    !< Print warning message
                ! end if
            type is (mesh_2D_Euler_heterog_c)
                !> Heterogeneous case computations for y-direction
                Courant_max=props_obj%flux_cent_y(1)*time_step/(&
                    props_obj%porosity(1)*mesh%Delta_y(1))  !< Initialize with first cell
                beta_max=props_obj%disp_cent(1)*time_step*mesh%sq_hypot(1)/(&
                    props_obj%porosity(1)*mesh%get_cell_size(1)**2)  !< Initialize Peclet with first cell value: Pe = Δx/α_L
                do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                    Courant=props_obj%flux_cent_y(i)*time_step/(&
                        props_obj%porosity(i)*mesh%Delta_y(i))  !< Compute local Courant
                    beta=props_obj%disp_cent(i)*time_step*mesh%sq_hypot(i)/(&
                        props_obj%porosity(i)*mesh%get_cell_size(i)**2)
                    if (Courant>Courant_max) then           !< Check if current exceeds maximum
                        Courant_max=Courant                 !< Update maximum
                    end if
                    if (beta>beta_max) then                 !< Check if current exceeds maximum
                        beta_max=beta                       !< Update maximum
                    end if
                end do
                this%beta=beta_max                     !< Store maximum beta number
                this%Courant_y=Courant_max                    !< Store maximum Courant number for y-direction
                if (this%Courant_y>1d0-this%Courant) then                  !< Check Courant condition
                    write(*,'(A,ES10.3,A,ES10.3,A)') ' WARNING: Courant_y = ', this%Courant_y, &
                        ' > 1 - Courant_x = ', 1d0-this%Courant, ' (Co_x + Co_y > 1)'
                    write(*,'(A)') '   Reduce Delta_t or refine mesh in y-direction'
                end if
                if (this%beta>0.5d0+1d-10) then               !< Check dispersion condition (with tolerance)
                    write(*,'(A,ES10.3,A)') ' WARNING: Dispersion number beta = ', this%beta, ' > 0.5'
                    write(*,'(A)') '   Reduce Delta_t or refine mesh'
                end if
                this%cell_Peclet=this%compute_cell_Peclet_2D(&
                    mesh%Delta_x(1),mesh%Delta_y(1),props_obj%long_dispersivity)  !< Initialize with first cell: Pe = Δx/α_L
                do i=2,mesh%Num_targets-mesh%targets_flag   !< Loop over all cells
                    Peclet=this%compute_cell_Peclet_2D(&
                        mesh%Delta_x(i),mesh%Delta_y(i),props_obj%long_dispersivity)
                    if (Peclet>this%cell_Peclet) then             !< Check if current exceeds maximum
                        this%cell_Peclet=Peclet                   !< Update maximum 
                    end if
                end do
                if (this%cell_Peclet>2d0) then                   !< Check Peclet condition
                    write(*,'(A,ES10.3,A)') ' WARNING: Cell Peclet number Pe = ', this%cell_Peclet, ' > 2.0'
                    write(*,'(A)') '   Refine mesh or increase dispersivity'
                end if
            end select
        end subroutine

        subroutine compute_cell_Peclet_1D(this,Delta_x,alpha_L)
            implicit none
            class(stab_params_tpt_1D_c) :: this                    !< Stability parameters object (modified)
            real(kind=8), intent(in) :: Delta_x                   !< Grid spacing Δx [L]
            real(kind=8), intent(in) :: alpha_L                    !< Longitudinal dispers            
            this%cell_Peclet=Delta_x/alpha_L  !< Peclet number: Pe = Δx/α_L
        end subroutine

        function compute_cell_Peclet_2D(this,Delta_x,Delta_y,alpha_L) result(cell_Peclet)
            implicit none
            class(stab_params_tpt_2D_c) :: this                    !< Stability parameters object (modified)
            real(kind=8), intent(in) :: Delta_x                   !< Grid spacing Δx [L]
            real(kind=8), intent(in) :: Delta_y                   !< Grid spacing Δy [L]
            !real(kind=8), intent(in) :: flux(:,:)                    !< Flux in x and y-direction at centers [L/T]    
            real(kind=8), intent(in) :: alpha_L                    !< Longitudinal dispersivity [L]            
            real(kind=8) :: cell_Peclet
            cell_Peclet=sqrt(Delta_x*Delta_y)/alpha_L  !< Peclet number: Pe = sqrt(Δx·Δy)/α_L
        end function
            ! real(kind=8) :: Peclet,Peclet_max                     !< Local and maximum Peclet numbers        
            ! Peclet_max=sqrt(flux(1,1)**2+flux(1,2)**2)*sqrt(Delta_x(1)*Delta_y(1))/&
            !     alpha_L  !< Peclet number in first cell: Pe = |q|·sqrt(Δx·Δy)/D
            ! do i=2,size(flux,1)  !< Loop over all cells
            !     Peclet=sqrt(flux(i,1)**2+flux(i,2)**2)*sqrt(Delta_x(i)*Delta_y(i))/&
            !         alpha_L  !< Peclet number in cell i: Pe = |q|·sqrt(Δx·Δy)/D
            !     if (Peclet>Peclet_max) then  !< Check if current exceeds maximum
            !         Peclet_max=Peclet  !< Update maximum
            !     end if
            ! end do
            !this%Peclet_y=Delta_y/alpha_L  !< Peclet number in y-direction: Pe = Δy/α_L
        ! end subroutine
end module