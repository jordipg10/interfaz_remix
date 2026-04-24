!> \file compute_trans_mat_tpt.f90
!> \brief Computes the transition matrix T for stationary advection-dispersion.
!> \details
!> Assembles the tridiagonal, negative semi-definite transition
!> matrix \f$ T \f$ from the finite-volume discretisation of the
!> 1D advection-dispersion equation \f$ 0 = T\,c + g \f$.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup discretization
!> \see transport_m, compute_trans_mat_diff

subroutine compute_trans_mat_tpt(this)                                          !< Compute transport matrix for advection-dispersion equation
!> Mathematical background: T: transition matrix (tridiagonal, negative semi-definite) !< Transport operator matrix T for spatial discretization
!> Mass conservation constraint: rows sum = 0 if r=0                       !< Row sums equal zero for conservative transport (no sources/sinks)
!> Governing equation: 0=T*c+g                                             !< Discretized ADE: T*concentration + source_term = 0
    use transport_m, only: transport_1D_c                                   !< Import 1D transport class
    use spatial_discr_1D_m                                                  !< Import spatial discretization modules
    use vectors_m                                                           !< Import vector operations
    implicit none                                                           !< Require explicit variable declaration
    
    class(transport_1D_c) :: this                                          !< Transport object containing spatial discretization and properties
    
    real(kind=8) :: sign_flux                                               !< Sign of flux (positive = downstream, negative = upstream)
    integer(kind=4) :: i,n                                                  !< Loop index and number of grid points
    
    n=this%spatial_discr%Num_targets                                        !< Get number of target points (computational nodes)

    if (this%spatial_discr%adapt_ref.eq.1) then                            !< Check if adaptive mesh refinement is active
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super) !< Deallocate existing matrix storage for mesh adaptation
        call this%allocate_trans_mat()                                     !< Reallocate transport matrix with new mesh dimensions
    end if                                                                  !< End adaptive mesh handling

    select type (mesh=>this%spatial_discr)                                 !< Polymorphic selection based on mesh type
    type is (mesh_1D_Euler_homog_c)                                        !< Homogeneous Eulerian mesh (uniform grid spacing)
        if (mesh%targets_flag.eq.0) then                                   !< Cell-centered discretization (unknowns at cell centers)
            if (mesh%scheme.eq.1) then                                     !< Traditional centered finite differences scheme
                this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2) + &
                    this%tpt_props_heterog%flux_cent(2:n)/(2*mesh%Delta_x)
                this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) - &
                    this%tpt_props_heterog%flux_cent(1:n-1)/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.2) then !> proposed by Petchamé & Carrera (2024)
                this%trans_mat%super(1)=this%tpt_props_heterog%disp_int(1)/(mesh%Delta_x**2) - &
                    this%tpt_props_heterog%flux_int(2)/(2*mesh%Delta_x)
                do i=2,n-1
                    this%trans_mat%sub(i-1)=this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) + &
                        this%tpt_props_heterog%flux_int(i)/(2*mesh%Delta_x)
                    this%trans_mat%super(i)=this%tpt_props_heterog%disp_int(i)/(mesh%Delta_x**2) - &
                        this%tpt_props_heterog%flux_int(i+1)/(2*mesh%Delta_x)
                end do
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%disp_int(n)/(mesh%Delta_x**2) + &
                    this%tpt_props_heterog%flux_int(n)/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.3) then !> upwind
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1))
                    this%trans_mat%sub=this%tpt_props_heterog%disp_cent(2:n)/(mesh%Delta_x**2) + &
                        ((sign_flux+1d0)/2)*this%tpt_props_heterog%flux_cent(2:n)/mesh%Delta_x
                    this%trans_mat%super=this%tpt_props_heterog%disp_cent(1:n-1)/(mesh%Delta_x**2) + &
                        ((sign_flux-1d0)/2)*this%tpt_props_heterog%flux_cent(1:n-1)/mesh%Delta_x
                end if
            else
                error stop "Scheme not implemented yet"
            end if
        end if
    end select
    this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term
    this%trans_mat%diag(2:n-1)=this%trans_mat%diag(2:n-1) - this%trans_mat%sub(1:n-2) - this%trans_mat%super(2:n-1)
end subroutine