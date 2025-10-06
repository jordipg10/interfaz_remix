!subroutine compute_char_params(this,k)
!>    use PDE_transient_m
!>    implicit none
!>    class(PDE_1D_transient_c) :: this
!>    integer(kind=4), intent(in), optional :: k
!>    !class(props_c), intent(in) :: props
!>    select type (props=>this%props)
!>    type is (tpt_props_homog_c)
!>        !this%char_params%disp_time=(this%spatial_discr%measure**2)/props%dispersion
!>        this%char_params%char_time=props%dispersion*props%porosity/(props%flux**2)
!>        this%char_params%char_measure=props%dispersion/props%flux
!>    end select
!end subroutine