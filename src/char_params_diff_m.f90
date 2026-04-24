!> \file char_params_diff_m.f90
!> \brief Diffusion characteristic parameters module.
!> \details
!> Extends `char_params_c` to compute characteristic time and length scales
!> for diffusion problems on radial meshes with homogeneous properties.
!>
!> \par Characteristic Scales (radial, homogeneous):
!> \f[
!>   L_c = r_{\max}, \quad t_c = \frac{\phi \, L_c^2}{D}
!> \f]
!> where \f$\phi\f$ is porosity, \f$D\f$ is diffusion coefficient, and
!> \f$r_{\max}\f$ is the outer radius.
!>
!> \see char_params_m, char_params_tpt_m, char_params_flow_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization
module char_params_diff_m
    use char_params_m, only: char_params_c
    use properties_m, only: props_c
    use spatial_discr_m, only: spatial_discr_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use diff_props_heterog_m, only: diff_props_heterog_1D_c
    use diffusion_m, only: diffusion_1D_c
    implicit none
    private
    !> \brief Diffusion characteristic parameters subclass.
    !> \details Computes \f$ t_c \f$ and \f$ L_c \f$ for diffusion problems.
    type, public, extends(char_params_c) :: char_params_diff_c
    contains
        procedure :: compute_char_params=>compute_char_params_diff  !< Compute diffusion characteristic scales
    end type

    contains

    !> \brief Compute characteristic parameters for diffusion.
    !> \details For homogeneous diffusion on a radial mesh:
    !> \f$ L_c = r_{\max} \f$ and \f$ t_c = \phi L_c^2 / D \f$.
    !> \param[in,out] this      Diffusion characteristic parameters object
    !> \param[in]     props_obj Physical properties (expects diff_props_heterog_1D_c)
    !> \param[in]     mesh      Spatial discretization (expects spatial_discr_rad_c)
    subroutine compute_char_params_diff(this,props_obj,mesh)
        class(char_params_diff_c) :: this
        class(props_c), intent(in) :: props_obj
        class(spatial_discr_c), intent(in) :: mesh
        select type (props_obj)
        type is (diff_props_heterog_1D_c) !> diffusion properties
            select type (mesh)
            type is (spatial_discr_rad_c) !> radial mesh
                if (props_obj%homog_flag.eqv..true.) then !> homogeneous properties
                    this%char_length=mesh%r_max !> r_c=r_ext
                    this%char_time=props_obj%porosity(1)*(this%char_length**2)/props_obj%diff_int(1) !> t_c=phi*r_c^2/D
                end if
            end select
        end select
    end subroutine

end module char_params_diff_m