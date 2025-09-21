module diff_stab_params_m
    use stability_parameters_m
    use diff_props_heterog_m
    use spatial_discr_rad_m
    use spatial_discr_1D_m
    use vectors_m
    implicit none
    save
    type, public, extends(stab_params_c) :: stab_params_diff_c !> 1D diffusion equation stability parameters subclass
        real(kind=8) :: beta !> dispersion stability parameter (beta=D*Delta_t/phi*Delta_x^2)
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_diff
    end type
    
    contains
        subroutine compute_stab_params_diff(this,props_obj,mesh,time_step)
            implicit none
            class(stab_params_diff_c) :: this
            class(props_c), intent(in) :: props_obj
            class(spatial_discr_c), intent(in) :: mesh
            !class(time_discr_c), intent(in) :: time_discr
            real(kind=8), intent(in) :: time_step !> time step
            
            integer(kind=4) :: i
            real(kind=8) :: beta_max,beta !> auxiliary variables for computing beta
                        
            select type (props=>props_obj)
            class is (diff_props_heterog_c)
                beta_max=props%dispersion(1)*time_step/(props%porosity(1)*mesh%get_cell_size(1)**2)
                do i=2,mesh%Num_targets-mesh%targets_flag
                    beta=props%dispersion(i)*time_step/(props%porosity(i)*mesh%get_cell_size(i)**2)
                    if (beta>beta_max) then
                        beta_max=beta
                    end if
                end do
                this%beta=beta_max
            end select
            if (this%beta>5d-1) then
                print *, "Unstable transport", this%beta
                !error stop "Dispersion condition violated"
            end if
        end subroutine
        
        
end module