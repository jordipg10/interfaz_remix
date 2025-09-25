module char_params_diff_m
    use char_params_m, only: char_params_c
    use properties_m, only: props_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c, spatial_discr_c
    use diff_props_heterog_m, only: diff_props_heterog_c
    use diffusion_m, only: diffusion_1D_c
    implicit none
    type, public, extends(char_params_c) :: char_params_diff_c !> diffusion characteristic parameters class
        !real(kind=8) :: char_time !> Characteristic time for diffusion
        !real(kind=8) :: char_measure !> Characteristic length for diffusion
    contains
        procedure, public :: compute_char_params=>compute_char_params_diff
    end type

    contains

    subroutine compute_char_params_diff(this,props_obj,mesh)
        class(char_params_diff_c) :: this
        class(props_c), intent(in) :: props_obj
        class(spatial_discr_c), intent(in) :: mesh
        select type (props_obj)
        type is (diff_props_heterog_c)
            select type (mesh)
            type is (spatial_discr_rad_c)
                if (props_obj%homog_flag.eqv..true.) then
                    this%char_length=mesh%r_max !> r_c=r_ext
                    this%char_time=props_obj%porosity(1)*(this%char_length**2)/props_obj%dispersion(1) !> t_c=phi*r_c^2/D
                end if
            end select
        end select
    end subroutine

end module char_params_diff_m