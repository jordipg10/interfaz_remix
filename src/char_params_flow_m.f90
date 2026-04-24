!> \file char_params_flow_m.f90
!> \brief Flow characteristic parameters module.
!> \details
!> Extends `char_params_c` to compute characteristic time, length, head, and
!> recharge scales for confined groundwater flow problems on radial meshes.
!>
!> \par Characteristic Scales (radial, homogeneous, confined):
!> \f[
!>   L_c = r_{\max}, \quad t_c = \frac{S \, L_c^2}{T}, \quad w_c = \frac{S \, h_c}{t_c}
!> \f]
!> where \f$S\f$ is storativity, \f$T\f$ is transmissivity, and \f$h_c\f$ is the
!> characteristic head chosen from boundary conditions.
!>
!> \see char_params_m, char_params_diff_m, char_params_tpt_m
!> \author Jordi
!> \date Unknown
!> \ingroup discretization
module char_params_flow_m
    use properties_m, only: props_c
    use char_params_m, only: char_params_c
    use flow_props_heterog_m, only: flow_props_heterog_c, flow_props_heterog_conf_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_m, only: spatial_discr_c
    use BCs_m, only: BCs_1D_c
    implicit none
    save
    !> \brief Flow characteristic parameters subclass.
    !> \details Adds characteristic head \f$ h_c \f$ and computes flow scaling.
    type, public, extends(char_params_c) :: char_params_flow_c
        real(kind=8) :: char_head  !< [L] Characteristic hydraulic head
    contains
        procedure :: compute_char_params=>compute_char_params_flow  !< Compute flow characteristic scales
        procedure :: set_char_head  !< Set characteristic head from boundary conditions
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
        !subroutine compute_char_params(this,props_obj,mesh)
        !    import char_params_c
        !    import props_c
        !    import spatial_discr_c
        !    implicit none
        !    class(char_params_c) :: this
        !    class(props_c), intent(in) :: props_obj
        !    class(spatial_discr_c), intent(in) :: mesh
        !end subroutine
    end interface
    
    contains
    
        subroutine compute_char_params_flow(this,props_obj,mesh)
            implicit none
            class(char_params_flow_c) :: this
            class(props_c), intent(in) :: props_obj
            class(spatial_discr_c), intent(in) :: mesh
            select type (props_obj)
            type is (flow_props_heterog_conf_c)
                select type (mesh)
                type is (spatial_discr_rad_c)
                    if (props_obj%homog_flag.eqv..true.) then
                        this%char_length=mesh%r_max !> r_c=r_ext
                        this%char_time=props_obj%storativity(1)*(this%char_length**2)/props_obj%transmissivity(1) !> t_c=S*r_c^2/T
                        this%char_rech=props_obj%storativity(1)*this%char_head/this%char_time !> w_c=S*h_c/t_c
                    end if
                end select
            end select
        end subroutine
        
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
        
        subroutine set_char_head(this,BCs)
            implicit none
            class(char_params_flow_c) :: this
            class(BCs_1D_c), intent(in) :: BCs !> characteristic head
            if (BCs%head_inf<=BCs%head_out) then
                this%char_head=BCs%head_out !> sets characteristic head to outflow head
            else
                this%char_head=BCs%head_inf !> sets characteristic head to inflow head
            end if
        end subroutine
        
        !subroutine set_char_head_rech(this,BCs,rech,
        !    implicit none
        !    class(char_params_flow_c) :: this
        !    class(BCs_1D_c), intent(in) :: BCs !> characteristic head
        !    if (BCs%head_inf<=BCs%head_out) then
        !        this%char_head=BCs%head_out !> sets characteristic head to outflow head
        !    else
        !        this%char_head=BCs%head_inf !> sets characteristic head to inflow head
        !    end if
        !end subroutine
end module