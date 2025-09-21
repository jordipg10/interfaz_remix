module stab_params_flow_m
    use stability_parameters_m, only: stab_params_c
    use properties_m, only: props_c
    use flow_props_heterog_m, only: flow_props_heterog_c, flow_props_heterog_conf_c
    use spatial_discr_m, only: spatial_discr_c
    implicit none
    save
    type, public, extends(stab_params_c) :: stab_params_flow_c !> 1D flow equation stability parameters subclass
        real(kind=8) :: beta !> flow stability parameter (beta=K*Delta_t/S_s*Delta_x^2 or beta=T*Delta_t/S*Delta_x^2)
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_flow
    end type
    
    contains
        subroutine compute_stab_params_flow(this,props_obj,mesh,time_step)
            implicit none
            class(stab_params_flow_c) :: this
            class(props_c), intent(in) :: props_obj
            !real(kind=8), intent(in) :: mesh_size
            class(spatial_discr_c), intent(in) :: mesh
            real(kind=8), intent(in) :: time_step
                        
            integer(kind=4) :: i
            real(kind=8) :: beta,beta_max
            
            select type (props=>props_obj)
            type is (flow_props_heterog_conf_c)
                this%Delta_t_crit=minval(props%storativity)*mesh%get_cell_size()**2/(2d0*maxval(props%transmissivity))
                beta_max=props%transmissivity(1)*time_step/(props%storativity(1)*mesh%get_cell_size(1)**2)
                do i=2,mesh%Num_targets-mesh%targets_flag
                    beta=props%transmissivity(i)*time_step/(props%storativity(i)*mesh%get_cell_size(i)**2)
                    if (beta>beta_max) beta_max=beta
                end do
            type is (flow_props_heterog_c)
                beta_max=props%hydr_cond(1)*time_step/(props%spec_stor(1)*mesh%get_cell_size(1)**2)
                do i=2,mesh%Num_targets-mesh%targets_flag
                    beta=props%hydr_cond(i)*time_step/(props%spec_stor(i)*mesh%get_cell_size(i)**2)
                    if (beta>this%beta) this%beta=beta
                end do
            end select
            this%beta=beta_max
            if (this%beta>5d-1) then
                print *, "Unstable flow", this%beta
                !error stop "Dispersion condition violated"
            end if
        end subroutine
        
        
end module