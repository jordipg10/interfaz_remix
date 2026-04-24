!> \file char_params_tpt_m.f90
!> \brief Transport characteristic parameters module.
!> \details
!> Extends `char_params_c` to compute characteristic time and length scales
!> for advection-dispersion transport problems on radial meshes with
!> homogeneous properties.
!>
!> \par Characteristic Scales (radial, homogeneous):
!> \f[
!>   L_c = r_{\max}, \quad t_c = \frac{\phi \, L_c^2}{D}
!> \f]
!> where \f$\phi\f$ is porosity, \f$D\f$ is the dispersion (centred) coefficient,
!> and \f$r_{\max}\f$ is the outer radius.
!>
!> \see char_params_m, char_params_diff_m, char_params_flow_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization
module char_params_tpt_m
    use char_params_m, only: char_params_c
    use properties_m, only: props_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_m, only: spatial_discr_c
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c
    use transport_m, only: transport_1D_c
    implicit none
    private
    !> \brief Transport characteristic parameters subclass.
    !> \details Computes \f$ t_c \f$ and \f$ L_c \f$ for transport problems.
    type, public, extends(char_params_c) :: char_params_tpt_c
    contains
        procedure :: compute_char_params=>compute_char_params_tpt  !< Compute transport characteristic scales
    end type

    contains

    !> \brief Compute characteristic parameters for transport.
    !> \details For homogeneous transport on a radial mesh:
    !> \f$ L_c = r_{\max} \f$ and \f$ t_c = \phi L_c^2 / D \f$.
    !> \param[in,out] this      Transport characteristic parameters object
    !> \param[in]     props_obj Physical properties (expects tpt_props_heterog_1D_c)
    !> \param[in]     mesh      Spatial discretization (expects spatial_discr_rad_c)
    subroutine compute_char_params_tpt(this,props_obj,mesh)
        class(char_params_tpt_c) :: this
        class(props_c), intent(in) :: props_obj
        class(spatial_discr_c), intent(in) :: mesh
        
        select type (props_obj)
        type is (tpt_props_heterog_1D_c)
            select type (mesh)
            type is (spatial_discr_rad_c)
                if (props_obj%homog_flag.eqv..true.) then !> homogeneous properties
                    this%char_length=mesh%r_max !> r_c=r_ext
                    this%char_time=props_obj%porosity(1)*(this%char_length**2)/props_obj%diff_cent(1) !> t_c=phi*r_c^2/D
                end if
            end select
        end select
    end subroutine

end module char_params_tpt_m