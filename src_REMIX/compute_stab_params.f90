!subroutine compute_stab_params(this)
!>    use PDE_transient_m
!>    implicit none
!>    class(PDE_1D_transient_c) :: this
!>    !class(props_c), intent(in) :: props
!>    select type (props=>this%props)
!>    class is (diff_props_homog_c)
!>        select type (mesh=>this%spatial_discr)
!>        type is (mesh_1D_Euler_homog_c)
!>            select type (time=>this%time_discr)
!>            type is (time_discr_homog_c)
!>                beta=props%dispersion*time%Delta_t/(props%porosity*mesh%Delta_x**2)
!>            end select
!>        end select
!>    end select
!>    if (this%beta>5d-1) error stop "Unstable diffusion"
!end subroutine