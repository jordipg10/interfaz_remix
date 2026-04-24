!> \file char_params_m.f90
!> \brief Abstract module for characteristic (dimensionless) parameters.
!> \details
!> Defines the abstract superclass `char_params_c` for computing characteristic
!> (reference) scales used in non-dimensionalization of PDE operators.
!> Subclasses compute specific characteristic time, length, and recharge scales
!> for diffusion, transport, or flow problems.
!>
!> \par Characteristic Scales:
!> - Characteristic time \f$ t_c \f$ [T]
!> - Characteristic length \f$ L_c \f$ [L]
!> - Characteristic recharge \f$ w_c \f$ [L/T]
!>
!> \see char_params_diff_m, char_params_tpt_m, char_params_flow_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization

!> \brief Abstract characteristic parameters module.
module char_params_m
    use properties_m, only: props_c
    use spatial_discr_m, only: spatial_discr_c
    implicit none
    save
    private
    !> \brief Abstract superclass for characteristic parameters.
    !> \details Stores reference scales for non-dimensionalization. Subclasses
    !> implement `compute_char_params` for specific PDE types.
    type, public, abstract :: char_params_c
        real(kind=8) :: char_time   !< [T] Characteristic time scale
        real(kind=8) :: char_length !< [L] Characteristic length scale
        real(kind=8) :: char_rech   !< [L/T] Characteristic recharge/source scale
    contains
        !> \brief Compute characteristic parameters from physical properties and mesh.
        procedure(compute_char_params), public, deferred :: compute_char_params
    end type
    
    !type, public, extends(char_params_c) :: char_params_diff_c !> diffusion characteristic parameters class
    !contains
    !    procedure :: compute_char_params=>compute_char_params_diff
    !end type
    !
    !type, public, extends(char_params_c) :: char_params_tpt_c !> transport characteristic parameters class
    !contains
    !    procedure :: compute_char_params=>compute_char_params_tpt
    !end type
    
    abstract interface
        subroutine compute_char_params(this,props_obj,mesh)
            import char_params_c
            import props_c
            import spatial_discr_c
            implicit none
            class(char_params_c) :: this
            class(props_c), intent(in) :: props_obj
            class(spatial_discr_c), intent(in) :: mesh
        end subroutine
    end interface
    
    contains
    
        !subroutine compute_char_params_diff(this,props_obj,mesh)
        !    implicit none
        !    class(char_params_diff_c) :: this
        !    class(props_c), intent(in) :: props_obj
        !    class(spatial_discr_c), intent(in) :: mesh
        !    select type (props_obj)
        !    type is (diff_props_heterog_1D_c)
        !        select type (mesh)
        !        type is (spatial_discr_rad_c)
        !            if (props_obj%homog_flag.eqv..true.) then
        !                this%char_measure=mesh%radius !> r_c=r_ext
        !                this%char_time=props_obj%porosity(1)*(this%char_measure**2)/props_obj%dispersion(1) !> t_c=phi*r_c^2/D
        !            end if
        !        end select
        !    end select
        !end subroutine
        !
        !subroutine compute_char_params_tpt(this,props_obj,mesh)
        !    implicit none
        !    class(char_params_tpt_c) :: this
        !    class(props_c), intent(in) :: props_obj
        !    class(spatial_discr_c), intent(in) :: mesh
        !    select type (props_obj)
        !    type is (tpt_props_heterog_1D_c)
        !        if (props_obj%homog_flag.eqv..true.) then
        !            this%char_time=props_obj%dispersion(1)*props_obj%porosity(1)/(props_obj%flux(1)**2) !> t_c=D*phi/q^2
        !            this%char_measure=props_obj%dispersion(1)/props_obj%flux(1) !> L_c=D/q
        !        end if
        !    end select
        !end subroutine
end module