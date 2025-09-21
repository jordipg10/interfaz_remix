module stability_parameters_m
    use time_discr_m
    use properties_m
    implicit none
    save
    type, public, abstract :: stab_params_c !> stability parameters superclass
        real(kind=8) :: Delta_t_crit !> critical time step
        real(kind=8) :: Delta_t_D_crit !> dimensionless critical time step
    contains
        procedure(compute_stab_params), public, deferred :: compute_stab_params
        procedure, public :: compute_Delta_t_D_crit
    end type
    
    abstract interface        
        subroutine compute_stab_params(this,props_obj,mesh,time_step)
            import stab_params_c
            import props_c
            import spatial_discr_c
            import time_discr_c
            implicit none
            class(stab_params_c) :: this
            class(props_c), intent(in) :: props_obj
            class(spatial_discr_c), intent(in) :: mesh
            !class(time_discr_c), intent(in) :: time_discr
            real(kind=8), intent(in) :: time_step !> time step
        end subroutine
    end interface
    
    contains
    
        subroutine compute_Delta_t_D_crit(this,t_c)
            implicit none
            class(stab_params_c) :: this
            !real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: t_c !> characteristic time
            
            this%Delta_t_D_crit=this%Delta_t_crit/t_c
        end subroutine
        
end module