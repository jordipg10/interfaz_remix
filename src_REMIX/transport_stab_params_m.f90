module transport_stab_params_m
    use diff_stab_params_m
    use transport_properties_heterog_m
    use spatial_discr_1D_m
    implicit none
    save
    type, public, extends(stab_params_diff_c) :: stab_params_tpt_c !> 1D transport stability parameters subclass
        real(kind=8) :: Courant !> advection stability parameter (=q*Delta_t/(phi*Delta_x))
        real(kind=8) :: Peclet !> Pe=Delta_x/alpha_L
    contains
        procedure, public :: compute_stab_params=>compute_stab_params_tpt
    end type
    
    contains
        subroutine compute_stab_params_tpt(this,props_obj,mesh,time_step)
        !> Computes transport stability parameters
        !! We assume that the mesh and time dicretisation are uniform
            implicit none
            class(stab_params_tpt_c) :: this !> stability parameters object
            class(props_c), intent(in) :: props_obj !> properties object
            class(spatial_discr_c), intent(in) :: mesh
            !class(time_discr_c), intent(in) :: time_discr
            real(kind=8), intent(in) :: time_step !> time step
            
            real(kind=8) :: D,phi,q !> dispersion coefficient, porosity, flux
            real(kind=8) :: Courant,Peclet,Courant_max,Peclet_max !> auxiliary variables
            real(kind=8), parameter :: epsilon=1d-12 !> tolerance
            integer(kind=4) :: i
            
            call compute_stab_params_diff(this,props_obj,mesh,time_step) !> computes diffusion stability parameters
            
            select type (props_obj)
            type is (tpt_props_heterog_c)
                call are_tpt_props_homog(props_obj)
                if (props_obj%homog_flag.eqv..true.) then !> homogeneous properties
                    phi=props_obj%porosity(1)
                    D=props_obj%dispersion(1)
                    q=props_obj%flux_cent(1)
                    this%Delta_t_crit=phi*mesh%get_cell_size()**2/(2d0*D)
                    this%Courant=abs(q)*time_step/(phi*mesh%get_cell_size())
                    if (this%Courant>1d0) then
                        print *, this%Courant
                        print *, "Courant condition violated"
                        !error stop  "Courant condition violated"
                    end if
                    this%Peclet=abs(q)*mesh%get_max_cell_size()/D
                    if (this%Peclet>2d0) then
                        print *, this%Peclet
                        print *, "Peclet condition violated"
                        !error stop  "Peclet condition violated"
                    end if
                else !> heterogeneous properties
                    !phi_max=maxval(props_obj%porosity)
                    !phi_min=minval(props_obj%porosity)
                    !!D_max=maxval(props_obj%dispersion)
                    !D_min=minval(props_obj%dispersion)
                    !q_max=maxval(props_obj%flux)
                    !!q_min=minval(props_obj%flux)
                    !this%Delta_t_crit=phi_max*mesh_size**2/(2d0*D_min)
                    !this%Courant=q_max*time_step/(phi_min*mesh_size)
                    Courant_max=props_obj%flux_cent(1)*time_step/(props_obj%porosity(1)*mesh%get_cell_size(1))
                    do i=2,mesh%Num_targets-mesh%targets_flag
                        Courant=props_obj%flux_cent(i)*time_step/(props_obj%porosity(i)*mesh%get_cell_size(i))
                        if (Courant>Courant_max) then
                            Courant_max=Courant
                        end if
                    end do
                    this%Courant=Courant_max
                    if (this%Courant>1d0) then
                        print *, this%Courant
                        print *, "Courant condition violated"
                    end if
                    Peclet_max=mesh%get_cell_size(1)/props_obj%dispersivity(1)
                    do i=2,mesh%Num_targets-mesh%targets_flag
                        Peclet=mesh%get_cell_size(i)/props_obj%dispersivity(i)
                        if (Peclet>Peclet_max) then
                            Peclet_max=Peclet
                        end if
                    end do
                    this%Peclet=Peclet_max
                    if (this%Peclet>2d0) then
                        print *, this%Peclet
                        print *, "Peclet condition violated"
                    end if
                end if
            end select            
        end subroutine
end module