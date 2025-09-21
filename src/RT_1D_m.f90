!> 1D reactive transport module (Lagrangian version)
module RT_1D_m
    use chemistry_m
    use transport_transient_m, only: transport_1D_transient_c, time_discr_homog_c, norm_mat_inf, tridiag_matrix_c
    use transport_m, only: transport_1D_c
    implicit none
    save
    type, public :: RT_1D_c !> 1D reactive transport superclass
        type(chemistry_c) :: chemistry                          !> chemistry object
    contains
        !procedure, public :: set_transport
        procedure, public :: set_chemistry
        procedure, public :: solve_RT_1D_ideal_lump_Lagr
        procedure, public :: solve_RT_1D_ideal_lump_Euler
        procedure, public :: write_RT_1D
        procedure, public :: write_python
        procedure, public :: compute_Delta_t_crit_RT
        procedure, public :: check_Delta_t_RT
        !procedure, public :: write_transport_data
        procedure, public :: read_time_discretisation
    end type
!***************************************************************************************************************************************************!
    type,public,extends(RT_1D_c) :: RT_1D_stat_c !> 1D stationary reactive transport subclass
        type(transport_1D_c) :: transport                       !> stationary transport object
    contains
        procedure, public :: set_transport_stat
    end type
!***************************************************************************************************************************************************!
    type,public,extends(RT_1D_c) :: RT_1D_transient_c !> 1D transient reactive transport subclass
        type(transport_1D_transient_c) :: transport             !> transient transport object
        real(kind=8) :: Delta_t_crit                            !> critical time step
        integer(kind=4) :: int_method_chem_reacts               !> integration method chemical reactions
    contains
        procedure, public :: set_int_method_chem_reacts
        procedure, public :: set_transport_trans
        procedure, public :: move_particles
        procedure, public :: introduce_particles
    end type
!***************************************************************************************************************************************************!
    interface
        
        subroutine solve_RT_1D_ideal_lump_Euler(this,root)
            import RT_1D_c
            implicit none
            class(RT_1D_c) :: this
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        subroutine solve_RT_1D_ideal_lump_Lagr(this,root)
            import RT_1D_c
            implicit none
            class(RT_1D_c) :: this
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        subroutine write_RT_1D(this,root,path_py)
            import RT_1D_c
            implicit none
            class(RT_1D_c), intent(in) :: this
            !integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
            character(len=*), intent(in), optional :: path_py
        end subroutine
        
        subroutine write_python(this,path)
            import RT_1D_c
            implicit none
            class(RT_1D_c), intent(in) :: this
            character(len=*), intent(in) :: path !> path output files
        end subroutine
       
        
        subroutine compute_Delta_t_crit_RT(this)
            import RT_1D_c
            implicit none
            class(RT_1D_c) :: this
        end subroutine
        
       subroutine read_time_discretisation(this,root)
            import RT_1D_c
            class(RT_1D_c) :: this
            !integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
        end subroutine
        
       subroutine read_transport_data(this,unit,file_tpt,mixing_ratios)
            import RT_1D_c
            class(RT_1D_c) :: this
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: file_tpt
        end subroutine
        
        !subroutine write_transport_data(this,unit)
        !    import RT_1D_c
        !    class(RT_1D_c) :: this
        !    integer(kind=4), intent(in) :: unit
        !end subroutine

    end interface
    
    contains
        subroutine set_transport_stat(this,transport_obj)
            implicit none
            class(RT_1D_stat_c) :: this
            class(transport_1D_c), intent(in) :: transport_obj
            this%transport=transport_obj
        end subroutine
        
        subroutine set_transport_trans(this,transport_obj)
            implicit none
            class(RT_1D_transient_c) :: this
            class(transport_1D_transient_c), intent(in) :: transport_obj
            this%transport=transport_obj
        end subroutine
        
        subroutine set_chemistry(this,chemistry_obj)
            implicit none
            class(RT_1D_c) :: this
            class(chemistry_c), intent(in) :: chemistry_obj
            this%chemistry=chemistry_obj
        end subroutine
        
      
        
        subroutine check_Delta_t_RT(this)
            implicit none
            class(RT_1D_c) :: this    
            
            real(kind=8), parameter :: eps=1d-16
            
            select type (this)
            class is (RT_1D_transient_c)
                if (this%transport%time_discr%int_method.eq.1) then
                    call this%compute_Delta_t_crit_RT()
                    if (abs(this%Delta_t_crit)<eps) then
                        continue
                    else
                        select type (time=>this%transport%time_discr)
                        type is (time_discr_homog_c)
                            if (time%Delta_t>this%Delta_t_crit) then
                                print *, this%Delta_t_crit
                                error stop "Delta_t is larger than Delta_t_crit"
                            end if
                        end select
                    end if
                end if
            end select
        end subroutine
        
        subroutine set_int_method_chem_reacts(this,int_method)
            implicit none
            class(RT_1D_transient_c) :: this
            integer(kind=4), intent(in) :: int_method
            if (int_method<1 .or. int_method>3) error stop "Integration method for chemical reactions not implemented yet"
            this%int_method_chem_reacts=int_method
        end subroutine

        subroutine move_particles(this,k)
        !> Move particles to the new coordinates based on the current transport spatial discretisation
        !> This is used in the transient case to update particle positions after each time step.
            implicit none
            class(RT_1D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in), optional :: k !> time step index
            
            integer(kind=4) :: target_ind !> target index associated to the displaced particle
            integer(kind=4) :: num_exit !> number of particles that left the domain
            real(kind=8), allocatable :: new_coords(:) !> new coordinates of the particles
            integer(kind=4) :: i !> loop index
            logical :: exit_flag !> flag to check if the particle has left the domain
            type(aqueous_chemistry_c), allocatable :: old_particles(:) !> old particles
            integer(kind=4), allocatable :: old_dom_indices(:) !> old particles indices
            integer(kind=4), allocatable :: old_bd_indices(:) !> old particles indices
            
            num_exit=0 !> initialise number of particles that left the domain
            
            allocate(new_coords(this%transport%spatial_discr%dim))
            old_particles=this%chemistry%target_waters !> get old particles
            old_dom_indices=this%chemistry%dom_tar_wat_indices !> get old domain indices
            old_bd_indices=this%chemistry%bd_waters_indices !> get old boudnary indices
            do i=1,this%chemistry%num_target_waters
                new_coords = this%chemistry%target_waters(i)%solid_chemistry%tar%coord + &
                    this%transport%tpt_props_heterog%flux_cent(&
                    this%chemistry%target_waters(i)%solid_chemistry%tar%id)*&
                    this%transport%time_discr%get_Delta_t(k)/&
                    this%transport%tpt_props_heterog%porosity(this%chemistry%target_waters(&
                    i)%solid_chemistry%tar%id)
                !print *, this%chemistry%target_waters(this%chemistry%dom_tar_wat_indices(i))%solid_chemistry%tar%id
                ! new_coords=this%chemistry%target_waters(this%chemistry%dom_tar_wat_indices(i))%solid_chemistry%tar%coord + &
                !     0.5*(this%transport%spatial_discr%get_cell_size(i)+this%transport%spatial_discr%get_cell_size(i+1))
                target_ind=this%transport%spatial_discr%get_target_ind(new_coords)
                if (target_ind>0) then
                    call this%chemistry%target_waters(i)%set_solid_chemistry(&
                        this%chemistry%target_solids(target_ind))
                else
                    call this%transport%spatial_discr%check_exit(new_coords,exit_flag)
                    if (exit_flag) then
                        print *, "Particle left the domain at time step ", k
                        print *, "Particle index: ", i
                        num_exit=num_exit+1 !> increment number of particles that left the domain
                    else
                        error stop "Particle moved to a non-existing target"
                    end if
                end if
            end do
            if (num_exit>0) then
                call this%chemistry%allocate_target_waters(this%chemistry%num_target_waters-num_exit) !> allocate new particles array
                call this%chemistry%allocate_dom_tar_wat_indices(this%chemistry%num_target_waters_dom-num_exit) !> allocate new indices array
                this%chemistry%dom_tar_wat_indices=old_dom_indices(1:size(old_dom_indices)-num_exit) !> set new indices
                !call this%chemistry%allocate_bd_wat_indices(this%chemistry%num_target_waters_dom-num_exit) !> allocate new indices array
                this%chemistry%bd_waters_indices(2)=this%chemistry%num_target_waters !> set new boundary water index
                this%chemistry%target_waters=old_particles(1:this%chemistry%num_target_waters) !> set new particles
                !this%chemistry%target_waters(this%chemistry%bd_waters_indices)=old_particles(old_bd_indices) !> set old particles to the chemistry object
                !this%chemistry%target_waters(this%chemistry%dom_tar_wat_indices)=old_particles(old_dom_indices(1:size(old_dom_indices)-num_exit)) !> set old particles to the chemistry object
            end if
        end subroutine

        subroutine introduce_particles(this,k)
        !> Introduce new particles at the first target based on the inflow conditions
        !> This is used in the transient case to introduce new particles at each time step.
            implicit none
            class(RT_1D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in), optional :: k !> time step index

            type(aqueous_chemistry_c), allocatable :: old_particles(:) !> old particles
            integer(kind=4), allocatable :: old_dom_indices(:) !> old particles indices
            integer(kind=4), allocatable :: old_bd_indices(:) !> old particles indices

            old_particles=this%chemistry%target_waters !> get old particles
            old_dom_indices=this%chemistry%dom_tar_wat_indices !> get old domain indices
            old_bd_indices=this%chemistry%bd_waters_indices !> get old boundary indices
            call this%chemistry%allocate_target_waters(this%chemistry%num_target_waters+1) !> set old particles to the chemistry object
            call this%chemistry%allocate_dom_tar_wat_indices(this%chemistry%num_target_waters_dom+1) !> allocate new indices array
            this%chemistry%bd_waters_indices(2)=this%chemistry%num_target_waters !> set new boundary water index
            this%chemistry%dom_tar_wat_indices(1:this%chemistry%num_target_waters_dom-1)=old_dom_indices !> set new indices
            this%chemistry%dom_tar_wat_indices(this%chemistry%num_target_waters_dom)=this%chemistry%num_target_waters-1 !> set the last domain target water index
            this%chemistry%target_waters(1)=this%chemistry%target_waters_init(1) !> set the first target water
            this%chemistry%target_waters(2:this%chemistry%num_target_waters)=old_particles !> set displaced particles
            call this%chemistry%target_waters(1)%set_volume(this%transport%BCs%caudal_inf*this%transport%time_discr%get_Delta_t(k)) !> set volume of the first target water
        end subroutine
end module RT_1D_m